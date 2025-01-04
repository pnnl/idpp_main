"""
    idpp/probability/analysis.py

    Dylan Ross (dylan.ross@pnnl.gov)

    module with KDTree subclasses for probability analysis
"""


from typing import Any, Union, Tuple, Dict, Set, List, Optional
import os
import errno
import pickle
from dataclasses import dataclass
import json

import numpy as np
import numpy.typing as npt
from sklearn.neighbors import KDTree
from scipy.sparse import coo_array, csr_array

from idpp.db.util import IdPPdb


# set up some type aliases
type PropertyTree = Union[MzTree, RtTree, CcsTree, Ms2Tree]
# map cmpd_id (int) to matching cmpd_ids (set(int))
type QueryResult = Dict[int, Set[int]]  
# map adduct_id (int) to similarity contributions from other adducts
# similarity contributions: data arrays for COO matrix
type Similarities = Tuple[
    npt.NDArray[np.int32], 
    npt.NDArray[np.int32], 
    npt.NDArray[np.float32]
]


class MzTree(KDTree):
    """ 
    a KDTree subclass for querying m/z values
    
    Attributes
    ----------
    mzs : ``numpy.ndarray(float)``
        input array of m/zs
    cmpd_ids : ``numpy.ndarray(int)``
        input array of cmpd_ids 
    """

    def __init__(self, 
                 mzs: npt.NDArray[np.float64], 
                 cmpd_ids: npt.NDArray[np.int32],
                 leaf_size: int = 128
                 ) :
        """
        create a new MzTree instance from an array of m/zs

        Parameters
        ----------
        mzs : ``numpy.ndarray(float)``
            input array of m/zs
        cmpd_ids : ``numpy.ndarray(int)``
            input array of cmpd_ids
        leaf_size : ``int``
            TODO
        """
        # store m/zs
        self.mzs = mzs
        self.cmpd_ids = cmpd_ids
        # call parent class __init__
        super().__init__(mzs.reshape(-1, 1), leaf_size=leaf_size)


    def _ppm_to_tol(self,
                    ppm: float
                    ) -> npt.NDArray[np.float64] :
        """
        (internal) compute an array of absolute tolerances from self.mzs using the specified ppm

        Parameters
        ----------
        ppm : ``float``
            search tolerance ppm
        
        Returns
        -------
        tolerances : ``numpy.ndarray(float)``
            tolerances for self.mzs computed from the specified ppm
        """
        return self.mzs * ppm / 1e6
    
    def query_all(self, 
                  ppm: float
                  ) -> QueryResult :
        """
        Search all of self.mzs using tolerance computed from specified ppm, 
        returns query result
        Parameters
        ----------
        ppm : ``float``
            search tolerance ppm

        Returns
        -------
        result : ``QueryResult``
            dict mapping cmpd_ids to sets of matching cmpd_ids
        """
        # TODO: Wrap the query_all_gen method to collect and return all results, instead of duplicating
        #       the logic in here.
        # this should actually be called tols because they are tolerances computed from specified ppm
        ppms = self._ppm_to_tol(ppm)  
        result = {}
        for idx in range(len(self.cmpd_ids)):
            mz = self.mzs[idx]
            cmpd_id = self.cmpd_ids[idx]
            tol = ppms[idx]
            mset = set([self.cmpd_ids[midx] for midx in self.query_radius([[mz,],], tol)][0])
            if (match_set := result.get(cmpd_id)) is not None:
                match_set |= mset
            else:
                result[cmpd_id] = mset
        return result
    
    def query_all_gen(self, 
                      ppm: float):
        """
        just like query_all() method but yields one search result at a time         
        """
        tols = self._ppm_to_tol(ppm)
        # Iterate through unique compound IDs. For each compound ID separately query all associated m/z
        # against the rest of the tree. Build up a set of matching compound IDs for each unique compound
        # ID. Yield one query result (i.e. unique compound ID and set of matching compound IDs) at a time.
        for qry_cmpd_id in np.unique(self.cmpd_ids):
            idx = np.where(self.cmpd_ids == qry_cmpd_id)[0]
            match_set = set()
            for qmz, qtol in zip(self.mzs[idx], tols[idx]):
                match_set |= set([self.cmpd_ids[midx] for midx in self.query_radius([[qmz,],], qtol)][0])
            yield qry_cmpd_id, match_set

    def save(self,
             dir : str,
             dataset_id: int
             ) -> None :
        """
        save this `MzTree` instance to file, load again using the `load_tree` function

        Parameters
        ----------
        dir : ``str``
            directory to save the tree instance into
        dataset_id : ``int``
            dataset identifier, used to generate file name
        """
        with open(os.path.join(dir, f"MzTree_dsid={dataset_id}.pkl"), "wb") as pf:
            pickle.dump((self, self.mzs, self.cmpd_ids), pf)

    def load_attrs(self,
                   *attrs
                   ) -> None :
        """ load the extra attributes that didnt get pickled automatically """
        mzs, cmpd_ids = attrs
        self.mzs = mzs
        self.cmpd_ids = cmpd_ids


class CcsTree(KDTree):
    """ 
    a KDTree subclass for querying CCS values
    
    Attributes
    ----------
    ccss : ``numpy.ndarray(float)``
        input array of CCS values
    cmpd_ids : ``numpy.ndarray(int)``
        input array of cmpd_ids 
    """

    def __init__(self, 
                 ccss: npt.NDArray[np.float64], 
                 cmpd_ids: npt.NDArray[np.int32],
                 leaf_size: int = 256
                 ) :
        """
        create a new MzTree instance from an array of m/zs

        Parameters
        ----------
        ccss : ``numpy.ndarray(float)``
            input array of CCS values
        cmpd_ids : ``numpy.ndarray(int)``
            input array with corresponding commpound IDs
        leaf_size : ``int``, default=64
            TODO
        """
        # store CCS array and compound IDs
        self.ccss = ccss
        self.cmpd_ids = cmpd_ids
        # call parent class __init__
        super().__init__(self.ccss.reshape(-1, 1), leaf_size=leaf_size)

    def _percent_to_tol_all(self,
                            percent : float
                            ) -> npt.NDArray[np.float64] :
        """
        (internal) compute an array of absolute tolerances from self.ccss using the 
        specified percent

        Parameters
        ----------
        percent : ``float``
            search tolerance percent
        
        Returns
        -------
        tolerances : ``numpy.ndarray(float)``
            tolerances for self.ccss computed from the specified percent
        """
        return (percent / 100.) * self.ccss
    
    def query_radius_single(self,
                            ccs: float,
                            percent: float
                            ) -> Set[int] :
        """
        Query a single CCS value using a specified radius and return the set of adduct IDs 
        within that radius

        Parameters
        ----------
        ccs : ``float``
            query CCS value
        percent : ``float``
            query tolerance (as a percent)

        Returns
        -------
        adduct_ids : ``set(int)``
            matching adduct IDs
        """
        return set(
            self.cmpd_ids[
                super().query_radius([[ccs]], (percent / 100.) * ccs)[0]
            ]
        )

    def query_radius(self, 
                     percent : float
                     ) -> npt.NDArray[Any] :
        """
        Search all of self.ccs_qry against the KDTree using a radius (tolerance computed from 
        specified percent),returns an array of all matching indices for each element in 
        self.ccs_qry

        .. note:: 

            This method is a thin wrapper around the ``KDTree.query_radius(...)`` method, and
            it returns the same array of arrays where each index in the first array contains
            an array of matching indices from the query. This is not so useful for my ultimate
            goal of coordinating the query results across multiple trees. So instead of using
            this method directly, the ``query_all(...)`` method should be used instead which 
            will take the output from this method and convert that into a dictionary that
            will incorporate cmpd_id info as well.

        Parameters
        ----------
        percent : ``float``
            search tolerance percent

        Returns
        -------
        TODO
        """
        return super().query_radius(self.ccss.reshape(-1, 1), self._percent_to_tol_all(percent))
    
    def query_all(self, 
                  percent : float
                  ) -> QueryResult :
        """
        Search all of self.ccs_qry using tolerance computed from specified percent
        returns query result

        Parameters
        ----------
        percent : ``float``
            search tolerance percent

        Returns
        -------
        result : ``QueryResult``
            dict mapping cmpd_ids to sets of matching cmpd_ids
        """
        return {
            self.cmpd_ids[idx]: set([self.cmpd_ids[midx] for midx in matches]) 
            for idx, matches in enumerate(self.query_radius(percent))
        }   
    
    def save(self,
             dir : str,
             dataset_id: int
             ) -> None :
        """
        save this `CcsTree` instance to file, load again using the `load_ccs_tree` function

        Parameters
        ----------
        dir : ``str``
            directory to save the tree instance into
        dataset_id : ``int``
            dataset identifier, used to generate file name
        """
        with open(os.path.join(dir, f"CcsTree_dsid={dataset_id}.pkl"), "wb") as pf:
            pickle.dump((self, self.ccss, self.cmpd_ids), pf)

    def load_attrs(self,
                   *attrs
                   ) -> None :
        """ load the extra attributes that didnt get pickled automatically """
        ccss, cmpd_ids = attrs
        self.ccss = ccss
        self.cmpd_ids = cmpd_ids
    

class RtTree(KDTree):
    """ 
    a KDTree subclass for querying RT values
    
    Attributes
    ----------
    rts : ``numpy.ndarray(float)``
        input array of RTs
    cmpd_ids : ``numpy.ndarray(int)``
        input array of cmpd_ids 
    """

    def __init__(self, 
                 rts: npt.NDArray[np.float64], 
                 cmpd_ids: npt.NDArray[np.int32],
                 leaf_size: int = 256
                 ) :
        """
        create a new RtTree instance from an array of m/zs

        Parameters
        ----------
        rts : ``numpy.ndarray(float)``
            input array of RTs
        cmpd_ids : ``numpy.ndarray(int)``
            input array of cmpd_ids
        leaf_size : ``int``
            TODO
        """
        # store RT array and compound ids
        self.rts = rts
        self.cmpd_ids = cmpd_ids
        # call parent class __init__
        super().__init__(self.rts.reshape(-1, 1), leaf_size=leaf_size)

    def query_radius_single(self,
                            rt: float,
                            tol: float
                            ) -> Set[int] :
        """
        Query a single RT value using a specified radius and return the set of adduct IDs 
        within that radius

        Parameters
        ----------
        rt : ``float``
            query RT value
        tol : ``float``
            query tolerance

        Returns
        -------
        adduct_ids : ``set(int)``
            matching adduct IDs
        """
        return set(
            self.cmpd_ids[
                super().query_radius([[rt]], tol)[0]
            ]
        )

    def query_radius(self, 
                     tol : float
                     ) -> npt.NDArray[Any] :
        """
        Search all of self.rts against the KDTree using a radius (tolerance in min.), 
        returns an array of all matching indices for each element in self.rts

        .. note:: 

            This method is a thin wrapper around the ``KDTree.query_radius(...)`` method, and
            it returns the same array of arrays where each index in the first array contains
            an array of matching indices from the query. This is not so useful for my ultimate
            goal of coordinating the query results across multiple trees. So instead of using
            this method directly, the ``query_all(...)`` method should be used instead which 
            will take the output from this method and convert that into a dictionary that
            will incorporate cmpd_id info as well.

        Parameters
        ----------
        tol : ``float``
            search tolerance (in min.)

        Returns
        -------
        TODO
        """
        return super().query_radius(self.rts.reshape(-1, 1), tol)
    
    def query_all(self, 
                  tol : float
                  ) -> QueryResult :
        """
        Search all of self.rts using tolerance in min.
        returns search result

        Parameters
        ----------
        tol : ``float``
            search tolerance (in min.)

        Returns
        -------
        result : ``QueryResult``
            dict mapping cmpd_ids to sets of matching cmpd_ids
        """
        return {
            self.cmpd_ids[idx]: set([self.cmpd_ids[midx] for midx in matches]) 
            for idx, matches in enumerate(self.query_radius(tol))
        }   
    
    def save(self,
             dir : str,
             dataset_id: int
             ) -> None :
        """
        save this `RtTree` instance to file, load again using the `load_tree` function

        Parameters
        ----------
        dir : ``str``
            directory to save the tree instance into
        dataset_id : ``int``
            dataset identifier, used to generate file name
        """
        with open(os.path.join(dir, f"RtTree_dsid={dataset_id}.pkl"), "wb") as pf:
            pickle.dump((self, self.rts, self.cmpd_ids), pf)

    def load_attrs(self,
                   *attrs
                   ) -> None :
        """ load the extra attributes that didnt get pickled automatically """
        rts, cmpd_ids = attrs
        self.rts = rts
        self.cmpd_ids = cmpd_ids


def _make_int_xlog2x_lookup():
    """ 
    make a lookup table with x*log2(x) values that can be indexed by the 
    sum of fragment intensities (in integer representation)
    """
    zero_to_two = np.linspace(0, 2, int(2e6 + 1), endpoint=True)
    return np.concatenate([[0], (zero_to_two[1:] * np.log2(zero_to_two[1:]))]).astype(np.float32)


class Ms2Tree:
    """ 
    a class for querying MS2 spectra with similar interface to KDTree
    
    Attributes
    ----------
    TODO
    cmpd_ids : ``numpy.ndarray(int)``
        input array of cmpd_ids 
    """

    def __init__(self, 
                 frag_imzs: npt.NDArray[np.int32], 
                 frag_iis: npt.NDArray[np.int32], 
                 adduct_ids: npt.NDArray[np.int32], 
                 cmpd_ids: npt.NDArray[np.int32]
                 ) :
        """
        create a new instance of Ms2Tree from array of ms2 fragments and associated cmpd_ids
        
        Parameters
        ----------
        frag_imzs : ``numpy.ndarray(int)``
            input array of framemt m/zs (in integer representation)
        frag_iis : ``numpy.ndarray(int)``
            input array of fragment abundances (in integer representation)
        cmpd_ids : ``numpy.ndarray(int)``
            input array of cmpd_ids 
        """
        # store the input arrays
        # they must all be pre-sorted by frag_imz
        sort_idx = np.argsort(frag_imzs)
        self.frag_imzs = frag_imzs[sort_idx]
        self.frag_iis = frag_iis[sort_idx]
        self.adduct_ids = adduct_ids[sort_idx]
        
        # NOTE: This dict comprehension relies on the assumption that a single adduct_id will never
        #       map to more than one cmpd_id. I think this assumption works based on how the schema
        #       is structured, but probably a good idea to make that clear here in case we end up getting
        #       weird results from this at some point down the road.
        # map adduct ids to compound ids
        self.adduct_to_cmpd_id = {
            adduct_id: cmpd_id 
            for adduct_id, cmpd_id in zip(self.adduct_ids, cmpd_ids[sort_idx])
        } 
        # store the length of the input arrays
        self.n = len(self.adduct_ids)
        # initialize the xlog2x lookup table
        self._xlog2x_lookup = _make_int_xlog2x_lookup()
        # precomputed similarities matrix is initially None
        # dict mapping adduct_id to dicts that map other adduct_ids to similarity scores
        self.similarities: Similarities = {}

    @property
    def _unique_adduct_ids(self):
        return self.adduct_to_cmpd_id.keys()

    def _similarity_contribution(self, 
                                 ii1: int, ii2: int
                                 ) -> np.float32 :
        """ compute a pair-wise similarity contribution """ 
        return self._xlog2x_lookup[ii1 + ii2] - self._xlog2x_lookup[ii1] - self._xlog2x_lookup[ii2]

    def _compute_group_similarity_contribution(self, 
                                               group_iis: Tuple[List[int], List[int]]
                                               ) -> None :
        """ 
        compute the similarity contributions between an accumulated group of fragments, 
        add to self._similarities 
        """
        # first compute coordinates and data for the group (as lists)
        grp_coord0: List[np.int32] = []
        grp_coord1: List[np.int32] = []
        grp_data: List[np.float32] = []
        # count members in the group
        i = 0
        for aid_A, ii_A in zip(*group_iis):
            for aid_B, ii_B in zip(*group_iis):
                # compute contribution
                sim: np.float32 = self._similarity_contribution(ii_A, ii_B)
                # sorted indices (low, high)
                i0, i1 = (aid_A, aid_B) if aid_A < aid_B else (aid_B, aid_A)
                # accumulate into group coords and data
                grp_coord0.append(i0)
                grp_coord1.append(i1)
                grp_data.append(sim)
                i += 1
        # convert to arrays and return
        # NOTE: https://stackoverflow.com/questions/30012362/faster-way-to-convert-list-of-objects-to-numpy-array
        return (
            np.fromiter(grp_data, np.float32, i),
            (np.fromiter(grp_coord0, np.int32, i), np.fromiter(grp_coord1, np.int32, i))
        )

    def precompute_similarities(self, 
                                imz_tol: int = 2000,
                                debug: bool = False,
                                ) -> None :
        """ 
        precompute a matrix of similarities between all spectra stored in this object 
        
        Parameters
        ----------
        imz_tol : ``int``, default=2000
            specify the tolerance (in Da, integer representation) for combining fragment mzs, the default
            is 2000 which corresponds to 20 mDa, the same that was used in the publication for flash
            entropy searchs
        """
        # 0.02000
        #    2000
        # determine the maximum index in the similarities matrix
        # the shape of the similarities matrix is the max index + 1
        max_sim_dim = max(self.adduct_ids) + 1
        sim_shape = (max_sim_dim, max_sim_dim)
        #print(f"{sim_shape=}")
        # iterate through the fragments arrays and compute similarity contributions from them
        # set initial conditions
        current_imz: np.int32 = self.frag_imzs[0]
        current_adduct_id: np.int32 = self.adduct_ids[0]
        current_ii: np.int32 = self.frag_iis[0]
        # NOTE: Every time the mz of the group gets jumped to a higher value, the upper limit mz for the 
        #       group becomes the current_imz + 2 * imz_tol. This means that the current mz at the time 
        #       the group is jumped is the lower limit, the effective mean group mz is 
        #       the current mz + tolerance, and the upper limit is that mean mz + the tolerance. 
        upper_imz: np.int32 = np.int32(current_imz + 2 * imz_tol)
        # group accumulators
        group_iis: Tuple[List[int], List[int]] = ([current_adduct_id], [current_ii])
        # index starts at 1 since the 0 elements have already been consumed while setting up the 
        # initial conditions
        idx: int = 1
        self.similarities = csr_array(sim_shape, dtype=np.float32)
        while idx < self.n:
            current_imz: np.int32 = self.frag_imzs[idx]
            current_adduct_id: np.int32 = self.adduct_ids[idx]
            current_ii: np.int32 = self.frag_iis[idx]
            if debug and idx % 1000 == 0:
                # if debug is true, update progress every so often
                print(f"\r    progress {100 * (idx + 1) / self.n:6.2f} %", end="    ")
            if current_imz > upper_imz:
                # compute similarity contributions for values in group accumulators
                # convert to CSR, add to similarities matrix
                self.similarities += coo_array(
                    self._compute_group_similarity_contribution(group_iis), 
                    shape=sim_shape
                ).tocsr()
                # clear out group accumulator
                group_iis = ([], [])
                # jump upper_imz to current_imz + 2 * imz_tol
                upper_imz = np.int32(current_imz + 2 * imz_tol)
            else:
                # add current imz, ii to group accumulator
                group_iis[0].append(current_adduct_id)
                group_iis[1].append(current_ii)
            # increment index
            idx += 1
        # if there are any values left in the group accumulators at the end of the loop, process them
        if len(group_iis) > 0:
            self.similarities += coo_array(
                self._compute_group_similarity_contribution(group_iis), 
                shape=sim_shape
            ).tocsr()
        # divide all of the similarity scores by 2 at the end (because the spectra normalized to sum to 1) 
        self.similarities *= 0.5
        if debug:
            print("\r    progress 100.00 %")

    def query_all(self, 
                  similarity_threshold : float
                  ) -> QueryResult :
        """
        Return query results for all spectra stored in this object using
        a specified spectral entropy similarity threshold.

        .. note::

            Must use ``precompute_similarities(...)`` method before performing queries

        Parameters
        ----------
        similarity_threshold : ``float``
            spectral entropy similarity threshold

        Returns
        -------
        result : ``QueryResult``
            dict mapping cmpd_ids to sets of matching cmpd_ids
        """
        # TODO: This needs to return a dict mapping COMPOUND IDs to COMPOUND IDs, not adduct IDs. The point
        #       of query all is to see how MS/MS spectra can narrow down compound annotations, so while the 
        #       actual queries are run at the level of adducts, the results should ultimately be rolled up
        #       up to the compound level.
        
        # return {
        #     self.cmpd_ids[idx]: set([self.cmpd_ids[midx] for midx in matches]) 
        #     for idx, matches in enumerate(self.query_radius(tol))
        # }
        if self.similarities == {}:
            msg = "self._similarities is not set, use precompute_similarities method before this method"
            raise RuntimeError(msg)
        result: QueryResult = {}
        # iterate through the adduct IDs
        for aid_A in self._unique_adduct_ids:
            for aid_B in self._unique_adduct_ids:
                # sorted indices (low, high) for similarities matrix
                i0, i1 = (aid_A, aid_B) if aid_A < aid_B else (aid_B, aid_A)
                if self.similarities[i0, i1] >= similarity_threshold:
                    if aid_A in result:
                        result[aid_A].add(aid_B)
                    else:
                        result[aid_A] = {aid_B}
        return result

    def save(self,
             dir : str,
             dataset_id: int
             ) -> None :
        """
        save this `Ms2Tree` instance to file, load again using the `load_tree` function

        Parameters
        ----------
        dir : ``str``
            directory to save the tree instance into
        dataset_id : ``int``
            dataset identifier, used to generate file name
        """
        with open(os.path.join(dir, f"Ms2Tree_dsid={dataset_id}.pkl"), "wb") as pf:
            pickle.dump((self,), pf)
    
    def load_attrs(self,
                   *attrs
                   ) -> None :
        """ load the extra attributes that didnt get pickled automatically """
        # NOTE: This method isn't actually necessary for this object because it is not a subclass
        #       of KDTree and therefore does not have the funky pickling behavior, but adding it
        #       so that it has the same loading behavior as the other tree classes.


def load_tree(tree_file: str
              ) -> PropertyTree :
    """
    load an already constructed `MzTree`, `CcsTree`, `Ms2Tree` or `RtTree` instance from file

    Parameters
    ----------
    file : ``str``
        file name for saved instance

    Returns
    -------
    tree : ``PropertyTree``
        loaded tree instance
    """
    # check that the file exists
    if not os.path.isfile(tree_file):
            raise FileNotFoundError(errno.ENOENT, 
                                    os.strerror(errno.ENOENT), 
                                    tree_file)
    with open(tree_file, "rb") as pf:
        tree, *attrs = pickle.load(pf)
    tree.load_attrs(*attrs)
    return tree


@dataclass
class DatasetQueries:
    """ Store a set of queries needed for selecting a complete dataset for identification probability analysis """
    mz_qry: Tuple[str, str]
    rt_qry: str
    ccs_qry: str
    ms2_qry: Tuple[str, str]

    def to_json(self) -> str :
        """ return a string with JSON representation of this object """
        return json.dumps({
            "mz_qry": self.mz_qry,
            "rt_qry": self.rt_qry,
            "ccs_qry": self.ccs_qry,
            "ms2_qry": self.ms2_qry
        }, indent=2)


def construct_mz_trees(db: IdPPdb,
                       queries: DatasetQueries
                       ) -> Tuple[MzTree, MzTree] :
    """
    Select a dataset for identification probability analysis using a set of input queries, then construct
    and return corresponding instance of `MzTree` for performing the analysis.

    Parameters
    ----------
    db : ``IdPPdb``
        interface for IdPP database
    queries : ``idpp.probability.analysis.DatasetQueries``
        Instance of DatasetQueries dataclass containing dataset selection querys

    Returns
    -------
    tree_pos : ``MzTree``
    tree_neg : ``MzTree``
        instances of`MzTree` constructed using the input queries, one for positive adducts the other for negative
    """
    #print("constructing property trees ...")
    mzts = [None, None]
    # construct the tree
    for i in range(2):
        cmpd_ids, props = [], []
        for cmpd_id, prop in db.cur.execute(queries.mz_qry[i]):
            cmpd_ids.append(cmpd_id)
            props.append(prop)
        mzts[i] = MzTree(np.array(props), np.array(cmpd_ids))
    #print(f"    MzTree ... ok")
    # unpack the individual trees
    # return the tree instance
    return tuple(mzts)


def construct_ccs_tree(db: IdPPdb,
                       queries: DatasetQueries
                       ) -> CcsTree :
    """
    Select a dataset for identification probability analysis using a set of input queries, then construct
    and return corresponding instance of `CcsTree` for performing the analysis.

    Parameters
    ----------
    db : ``IdPPdb``
        interface for IdPP database
    queries : ``idpp.probability.analysis.DatasetQueries``
        Instance of DatasetQueries dataclass containing dataset selection querys

    Returns
    -------
    tree : ``CcsTree``
        instance of`CcsTree` constructed using the input query
    """
    print("constructing property trees ...")
    # construct the tree
    cmpd_ids, props = [], []
    for cmpd_id, prop in db.cur.execute(queries.ccs_qry):
        cmpd_ids.append(cmpd_id)
        props.append(prop)
    ccst = CcsTree(np.array(props), np.array(cmpd_ids))
    print(f"    CcsTree ... ok")
    # unpack the individual trees
    # return the tree instance
    return ccst


def construct_rt_tree(db: IdPPdb,
                      queries: DatasetQueries
                      ) -> RtTree :
    """
    Select a dataset for identification probability analysis using a set of input queries, then construct
    and return corresponding instance of `RtTree` for performing the analysis.

    Parameters
    ----------
    db : ``IdPPdb``
        interface for IdPP database
    queries : ``idpp.probability.analysis.DatasetQueries``
        Instance of DatasetQueries dataclass containing dataset selection querys

    Returns
    -------
    tree : ``RtTree``
        instance of`RtTree` constructed using the input query
    """
    print("constructing property trees ...")
    # construct the tree
    cmpd_ids, props = [], []
    for cmpd_id, prop in db.cur.execute(queries.rt_qry):
        cmpd_ids.append(cmpd_id)
        props.append(prop)
    rtt = RtTree(np.array(props), np.array(cmpd_ids))
    print(f"    RtTree ... ok")
    # unpack the individual trees
    # return the tree instance
    return rtt


def construct_ms2_tree(db: IdPPdb,
                       queries: DatasetQueries
                       ) -> Optional[Ms2Tree] :
    """
    Select a dataset for identification probability analysis using a set of input queries, then construct
    and return corresponding instance of `Ms2Tree` for performing the analysis.

    Parameters
    ----------
    db : ``IdPPdb``
        interface for IdPP database
    queries : ``idpp.probability.analysis.DatasetQueries``
        Instance of DatasetQueries dataclass containing dataset selection querys

    Returns
    -------
    tree : ``Ms2Tree`` or None
        instance of`Ms2Tree` constructed using the input query, or None if no spectra
    """
    ms2_qry_A, ms2_qry_B = queries.ms2_qry
    aid_to_cnt = {
        aid: cnt
        for aid, cnt in db.cur.execute(ms2_qry_A)
    }
    if sum(aid_to_cnt.values()) > 0:
        cmpd_ids, adduct_ids, frag_imzs, frag_iis = [], [], [], []
        res_B = db.cur.execute(ms2_qry_B).fetchall()
        for i, (cmpd_id, adduct_id, frag_imz, frag_ii) in enumerate(res_B):
            cmpd_ids.append(cmpd_id)
            adduct_ids.append(adduct_id)
            frag_imzs.append(frag_imz)
            frag_iis.append(frag_ii // aid_to_cnt[adduct_id])
        # the lists are all type int, convert to np.int32
        return Ms2Tree(
            np.fromiter(frag_imzs, np.int32, i + 1), 
            np.fromiter(frag_iis, np.int32, i + 1), 
            np.fromiter(adduct_ids, np.int32, i + 1), 
            np.fromiter(cmpd_ids, np.int32, i + 1)
        )
    else:
        return None


def construct_ms2_tree_for_adduct_ids(db: IdPPdb, 
                                      adduct_ids: int,
                                      precompute_similarities: bool = True
                                      ) -> Optional[Ms2Tree] :
    """
    Alternate function for constructing an Ms2Tree from spectra based on a specified list of adduct IDs
    If there are no MS/MS spectra associated with the adduct_ids, returns None

    Parameters
    ----------
    db : ``IdPPdb``
        interface for IdPP database
    adduct_ids : ``list(int)``
        adduct identifiers to attempt to gather spectra from 
    precompute_similarities : ``bool``, default=True
        flag specifying whether to precompute spectra similarity matrix
    
    Returns
    -------
    ms2t : ``Ms2Tree`` or ``None``
        instance of Ms2Tree if there were MS/MS spectra associated with the input adduct_id, or else None 
    """
    ms2_qry_A = """--sqlite3
    SELECT 
        adduct_id,
        COUNT(*) AS cnt 
    FROM 
        MS2Spectra 
        JOIN 
            Adducts USING(adduct_id)
    WHERE
        adduct_id IN ({})
    GROUP BY 
        adduct_id
    ;"""
    ms2_qry_B = """--sqlite3
    SELECT 
        cmpd_id, 
        adduct_id,
        frag_imz, 
        SUM(frag_ii)
    FROM 
        MS2Spectra 
        JOIN 
            Adducts USING(adduct_id)
        JOIN
            Compounds USING(cmpd_id)
        JOIN
            MS2Fragments USING(ms2_id)
    WHERE
        adduct_id IN ({})
    GROUP BY 
        adduct_id, 
        frag_imz
    ;"""
    aid_s = ",".join([str(_) for _ in adduct_ids])
    ms2t = construct_ms2_tree(
        db,
        DatasetQueries(
            mz_qry="none", 
            rt_qry="none",
            ccs_qry="none",
            ms2_qry=(ms2_qry_A.format(aid_s), ms2_qry_B.format(aid_s))
        )
    )
    if precompute_similarities and ms2t is not None:
        ms2t.precompute_similarities()
    return ms2t


def construct_property_trees(db: IdPPdb,
                             queries: DatasetQueries
                             ) -> Tuple[int, Tuple[MzTree, RtTree, CcsTree, Ms2Tree]] :
    """
    Select a dataset for identification probability analysis using a set of input queries, then construct
    and return corresponding instances of `MzTree`, `RtTree`, `CcsTree` and `Ms2Tree` for performing
    the analysis.

    Parameters
    ----------
    db : ``IdPPdb``
        interface for IdPP database
    queries : ``idpp.probability.analysis.DatasetQueries``
        Instance of DatasetQueries dataclass containing dataset selection querys

    Returns
    -------
    trees : ``tuple(MzTree, RtTree, CcsTree, Ms2Tree)``
        instances of `MzTree`, `RtTree`, `CcsTree` and `Ms2Tree` constructed from data fetched
        using the input query
    """
    # TODO: Remove most logic from this function and just use the individual functions above
    #       (construct_X_tree). No need to duplicate the logic.
    print("constructing property trees ...")
    # construct the m/z, RT, and CCS trees 
    # the logic is about the same for all of them
    mrc_trees = []
    for qry, t, lbl in zip([queries.mz_qry, queries.rt_qry, queries.ccs_qry], 
                           [MzTree, RtTree, CcsTree],
                           ["MzTree", "RtTree", "CcsTree"]):
        cmpd_ids, props = [], []
        for cmpd_id, prop in db.cur.execute(qry):
            cmpd_ids.append(cmpd_id)
            props.append(prop)
        mrc_trees.append(t(np.array(props), np.array(cmpd_ids)))
        print(f"    {lbl} ... ok")
    # unpack the individual trees
    mzt, rtt, ccst = mrc_trees
    # construct the MS2 tree, this works a bit differently than the other properties
    ms2_qry_A, ms2_qry_B = queries.ms2_qry
    aid_to_cnt = {
        aid: cnt
        for aid, cnt in db.cur.execute(ms2_qry_A)
    }
    cmpd_ids, adduct_ids, frag_imzs, frag_iis = [], [], [], []
    res_B = db.cur.execute(ms2_qry_B).fetchall()
    n = len(res_B)
    for i, (cmpd_id, adduct_id, frag_imz, frag_ii) in enumerate(res_B):
        cmpd_ids.append(cmpd_id)
        adduct_ids.append(adduct_id)
        frag_imzs.append(frag_imz)
        frag_iis.append(frag_ii // aid_to_cnt[adduct_id])
        print(f"\r    Ms2Tree: {100 * (i + 1) / n:6.2f} % ", end="")
        # # ! EARLY STOP FOR TESTING !
        # if i >= 100_000:
        #     break
    print("\n... ok")
    # the lists are all type int, convert to np.int32
    ms2t = Ms2Tree(
        np.fromiter(frag_imzs, np.int32, i + 1), 
        np.fromiter(frag_iis, np.int32, i + 1), 
        np.fromiter(adduct_ids, np.int32, i + 1), 
        np.fromiter(cmpd_ids, np.int32, i + 1)
    )
    # return the tree instances
    return mzt, rtt, ccst, ms2t

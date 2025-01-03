"""
    idpp/probability/analysis.py

    Dylan Ross (dylan.ross@pnnl.gov)

    module for probability analysis functions
"""


from typing import Optional

import numpy as np

from idpp.db.util import IdPPdb
from idpp.probability.trees import MzTree, RtTree, CcsTree, Ms2Tree, load_tree


def aggregate_query_results(db: IdPPdb,
                            ds_id: int, 
                            mzt: MzTree, 
                            mz_ppm: float,
                            rtt: Optional[RtTree] = None, 
                            rt_tol: Optional[float] = None,
                            ccst: Optional[CcsTree] = None, 
                            ccs_percent: Optional[float] = None,
                            ms2t: Optional[Ms2Tree] = None,
                            ms2_sim: Optional[float] = None,
                            ) -> None :
    """  
    Perform queries of the property trees using a set of defined search tolerances, 
    and store the aggregated identification counts and query information in the 
    IdPPdb database
    
    Parameters
    ----------
    db : ``IdPPdb``
        interface for interacting with the IdPP database
    ds_id : ``int``
        dataset identifier to specify which probability analysis dataset these 
        results correspond to
    mzt : ``MzTree``
        instance of MzTree for running queries based on m/z
    mz_ppm : ``float``
        m/z search tolerance in ppm
    rtt : ``RtTree``, optional
        instance of RtTree for running queries based on RT, do not consider RT if None
    rt_tol : ``float``, optional
        RT search tolerance, must be provided if ``rtt`` is not None
    ccst : ``CcsTree``, optional
        instance of CcsTree for running queries based on CCS, do not consider CCS if None
    ccs_percent : ``float``, optional
        CCS search tolerance in percent, must be provided if ``ccst`` is not None
    ms2t : ``Ms2Tree``, optional
        instance of Ms2Tree for running queries based on MS2, do not consider MS2 if None
    ms2_sim : ``float``, optional
        similarity threshold for querying MS2, must be provided if ``ms2t`` is not None
    """
    # validate: make sure if any of the optional tree instances are not None, that their corresponding 
    # search tolerance is also not None (for instance if rtt is not None then rt_tol must also not be None)
    # gather query results from property trees
    q_res = mzt.query_all(mz_ppm)
    # aggregate between query results from different trees
    # only aggregate compound ids that are present in all query results. 
    # the aggregated results consist of a numpy array of integers with length N
    # where N is the number of compound ids that are common to results from all 
    # the properties being considered (see examples below)
    # each element in the aggregated results array corresponds to the count of matched compound ids 
    # (from all properties) for each common compound id
    """
    # EXAMPLE
    # A query result maps query compound IDs (key) to sets of matching compound IDs (value)
    # query_result = {query_cmpd_id_1: {matched_cmpd_id_1, matched_cmpd_id_2, ...}, query_cmpd_id_2: {...}, ...}
    
    # query result from MzTree.query_all()
    q_result_mz = {
        69: {69, 420}
        123: {123, 456, 789, 1012},
        234: {234, 567, 890},
        555: {555, 666, 777}
    }

    # query result from RtTree.query_all()
    q_result_rt = {
        123: {123, 456, 788},
        234: {234, 566, 890, 999},
        420: {420, 69},
        444: {444, 545, 646}
    }

    # common compound IDs between the two result sets
    common_cmpd_ids = set(q_result_mz.keys()) & set(q_result_rt.keys())

    # -> in this example should be [123, 234], these are the only common keys

    # iterate through common (query) compound ids, count up the matches across specified properties
    counts = []
    for q_cmpd_id in common_cmpd_ids:
        counts.append(len(q_result_mz[q_cmpd_id] & q_result_rt[q_cmpd_id]))

    # -> in this example counts should be [2, 2]
    """
    
    # store in IdPPdb 
    # counts should be a list of ints after the results have been aggregated
    counts = [0, 1, 2]
    # convert counts list to numpy array, then use .tobytes() method on the array
    # this converts the array into a binary format that we can store directly in 
    # the database (this is done inline within the call to insert_analysis_result)
    db.insert_analysis_result(ds_id, np.array(counts).tobytes(), mz_ppm, 
                              rt_tol=rt_tol, ccs_tol=ccs_percent, ms2_tol=ms2_sim)


def _():


    db = IdPPdb(...)

    mzt = load_tree("path/to/mztree.pkl")
    rtt = load_tree("path/to/rttree.pkl")
    ccst = load_tree("path/to/ccstree.pkl")
    ms2t = load_tree("path/to/ms2tree.pkl")
    
    ds_id = 1


    for ppm in [1, 3, 5, 10, 30]:
        aggregate(db, ds_id, mzt, rtt, ccst, ms2t, ppm)

    for ppm in [1, 3, 5, 10, 30]:
        for rttol in [0.1, 0.5, 1.0]:
            aggregate(db, ds_id, mzt, rtt, ccst, ms2t, ppm, rt_tol=rttol)



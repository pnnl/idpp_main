"""
    idpp/subsetting/rdkfp_and_clust.py

    Dylan Ross (dylan.ross@pnnl.gov)

    functions for subsetting a database using a workflow that includes
    computing molecular fingerprints using RDKit, then dimensionality
    reduction and clustering. 
"""


from rdkit import Chem
from rdkit.Chem.rdmolops import RDKFingerprint
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import euclidean_distances


def fingerprint_and_cluster(smis, rdkfp_size, n_pca_comp, n_clust, 
                            random_state=420, include_pca_metadata=False):
    """
    Takes an array of input molecules (as SMILES structures), computes RDKit fingerprints,
    performs PCA on the fingerprint data for dimensionality reduction, then clusters the 
    PCA data (K-Means). Ultimately returns the cluster IDs for each molecule and the indices 
    of the central structures from each cluster.
    
    Parameters
    ----------
    smis : ``numpy.ndarray(float)``
        input molecules as array of SMILES structures, shape: (n_mol,)
    rdkfp_size : ``int``
        size of the RDKit fingerprint in bits (default size in the original RDKit function is 2048)
    n_pca_comp : ``int``
        number of PCA components to use for dimensionality reduction of fingerprint data
    n_clust : ``int``
        number of clusters to use for K-Means clustering on the PCA data
    random_state : ``int``, default=420
        pRNG seed for deterministic results in the PCA and K-Means steps
    include_pca_metadata : ``bool``, default=False
        if True, returns a dictionary of metadata from the PCA step including pca projections
        of the fingerprint data and explained variance ratios for all components, if False 
        return `None`. This data can be pretty large for big databases so consider omitting
        if storage/memory is a concern

    Returns
    -------
    cluster_ids : ``numpy.ndarray(int)``
        array of cluster IDs for all input molecules, shape: (n_mol,)
    cluster_center_idxs : ``numpy.ndarray(int)``
        array of indices of input molecules which represent the most central structure within
        each cluster, shape: (n_clust,) 
    pca_metadata : ``dict(...) or None``
        if `include_pca_metadata` param was True, return a dictionary with PCA projections (key='proj')
        for all compounds and explained variance ratios from all components (key='evrs'), otherwise return
        `None`
    """
    # compute fingerprints from SMILES structures, construct unpacked array of `char`s from bitstrings
    fps = np.array([[int(_) for _ in RDKFingerprint(Chem.MolFromSmiles(smi), fpSize=rdkfp_size)] for smi in smis], dtype=np.uint8)
    # fit the PCA, compute projections
    pca = PCA(n_components=n_pca_comp, svd_solver='full', random_state=random_state).fit(fps)
    fps_pca = pca.transform(fps)
    # run K-Means clustering on the fingerprint data
    kmc = KMeans(n_clusters=n_clust, random_state=random_state, n_init='auto').fit(fps_pca)
    # find the central structures for each cluster
    # compute distances of all molecules against all cluster centers
    # this will have shape (n_mol, n_clust)
    # get the index of the minimum distance within each column, these are the central molecules
    cluster_center_idxs = np.argmin(euclidean_distances(fps_pca, kmc.cluster_centers_), axis=0)
    # gather the PCA metadata if specified
    pca_metadata = {'proj': fps_pca, 'evrs': pca.explained_variance_ratio_} if include_pca_metadata else None
    # return values
    return kmc.labels_, cluster_center_idxs, pca_metadata

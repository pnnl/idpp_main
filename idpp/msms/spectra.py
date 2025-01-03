"""
    idpp/msms/spectra.py

    Dylan Ross (dylan.ross@pnnl.gov)

    module with functions for comparing or combining spectra
"""


from typing import List

import numpy as np
from numpy import typing as npt


def spec_norm(spectrum):
    """
    normalize spectrum such that all intensities sum to 1
    
    Parameters
    ----------
    spectrum : ``list(list(float))``
        MS/MS spectrum as list of [[m/z values...], [intensities...]]

    Returns
    -------
    spectrum : ``numpy.ndarray()``
        normalized MS/MS spectrum with same shape as input
    """
    mz, i = np.array(spectrum)
    return np.array([mz, i / sum(i)])


def spec_entropy(spectrum):
    """
    compute spectral entropy for single spectrum 
    
    Parameters
    ----------
    spectrum : ``list(list(float))``
        MS/MS spectrum as list of [[m/z values...], [intensities...]]
    
    Returns
    -------
    entropy : ``float``
        spectral entropy
    """
    mz, i = spectrum
    return -np.sum(i * np.log(i + 0.001))
    

def spec_combine(spectra: List[npt.NDArray],
                 weights: List[float],
                 mztol: float = 0.05
                 ) -> npt.NDArray :
    """ 
    combine multiple spectra into a single spectrum 
    
    Parameters
    ----------
    spectra : ``list(numpy.ndarray())``
        list of MS/MS spectra (2D arrays with shape (2, n_points)) to combine
    weights : ``list(float)``
        weights for each spectrum
    mztol : ``float``, default=0.05
        m/z tolerance for combining m/z peaks from different spectra

    Returns
    -------
    comb_spectrum : ``numpy.ndarray()``
        combined MS/MS spectrum (shape: (2, n_points))
    """
    assert len(spectra) == 2, "can only combine 2 spectra at a time"
    assert len(weights) == 2, "need exactly 2 weights"
    # unpack spectra
    (comb_mz, comb_i), (other_mz, other_i) = spectra
    comb_mz, comb_i = comb_mz.tolist(), comb_i.tolist()
    # apply weights
    comb_i = [weights[0] * _ for _ in comb_i] if weights[0] != 1. else comb_i
    other_i = [weights[1] * _ for _ in other_i] if weights[1] != 1. else other_i
    # iterate through both lists and combine if needed
    not_comb_mz, not_comb_i = [], [] 
    comb_n, other_n = len(comb_mz), len(other_mz)
    comb_idx, other_idx = 0, 0
    while comb_idx < comb_n and other_idx < other_n:
        omz, oi = other_mz[other_idx], other_i[other_idx]
        dmz = omz - comb_mz[comb_idx]
        if abs(dmz) <= mztol:
            comb_i[comb_idx] += oi
            comb_idx += 1
            other_idx += 1
        elif dmz < 0:
            not_comb_mz.append(omz)
            not_comb_i.append(oi)
            other_idx += 1
        else:
            comb_idx += 1
    # if there are still elements in other, add them
    while other_idx < other_n:
        comb_mz.append(other_mz[other_idx])
        comb_i.append(other_i[other_idx])
        other_idx += 1
    comb_mz += not_comb_mz
    comb_i += not_comb_i
    # sort by m/z
    idx = np.argsort(comb_mz)
    # return spectrum sorted by m/z
    # do not normalize
    return np.array(comb_mz)[idx], np.array(comb_i)[idx]


def spec_entropy_similarity(spectrum_A, spectrum_B):
    """ 
    pairwise spectral entropy based similarity as defined in the paper:
        https://www.nature.com/articles/s41592-021-01331-z
    
    Parameters
    ----------
    spectrum_A, spectrum_B : ``numpy.ndarray()``
        input MS/MS spectra (2D arrays with shape (2, n_points)) to compare
    
    Returns
    -------
    similarity : ``float``
        spectral entropy similarity score
    """
    spectrum_AB = spec_combine([spectrum_A, spectrum_B])
    s_AB, s_A, s_B = spec_entropy(spectrum_AB), spec_entropy(spectrum_A), spec_entropy(spectrum_B)
    return 1. - ((2. * s_AB - s_A - s_B) / np.log(4.))


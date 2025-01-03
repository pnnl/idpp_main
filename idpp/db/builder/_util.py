"""
    idpp/db/builder/_util.py

    Dylan Ross (dylan.ross@pnnl.gov)

    (internal) General utilities for building the database
"""


from typing import Tuple, Optional
import re

import numpy as np
from numpy import typing as npt


def parse_ce(ce_str: Optional[str]
             ) -> Optional[int] :
    """
    convert the various collision energy descriptions to an integer value where 
    possible otherwise return None 

    Parameters
    ----------
    ce_str : ``str`` or ``None``
        collision energy description, can be None

    Returns
    -------
    ce : ``int`` or ``None``
        if description was not None and was parsable, 
        then return the collision energy as an int
    """
    if ce_str is None:
        return None
    # stuff like "35 eV" or "40 V" or "22eV" or "54V"
    # possibly with decimal 
    # or even just a plain number
    # but no %
    # maybe with CE in the front
    # or NCE/(NCE) after 
    # wow would you look at that pattern lol
    m = re.fullmatch('(N*CE)*[- ]*([0-9]+([.][0-9]*)*)( *e*V)*( *[(]*NCE[)]*)*', ce_str)
    if m:
        # TODO: I have spent too long on the pattern trying to exclude % and can't figure 
        #       it out without breaking either the expected parsable or expected unparsable
        #       cases so I am adding an explicit check here to filter out things like 
        #       "50%" or "50 %". Maybe someone better at regex can bake that into the pattern
        #       and make this cleaner? The idpp.test.db.builder.mona.Test_ParseCe test case
        #       class has the expected parsable and unparsable inputs so that can be used 
        #       to check the pattern.
        if "%" in ce_str:
            return None
        # round to nearest int
        return int(round(float(m.group(2)), 0))    
    # didn't match anything we expected
    return None


def str_to_ms2(s: str
               ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]] :
    """
    converts flat str representation (space-separated peaks in format "{mz}:{intensity}") 
    to arrays of m/z and intensities

    Parameters
    ----------
    s : ``SpecStr``
        string form of spectrum

    Returns
    -------
    spectrum : ``Spec``
        spectrum as 2D array with m/z and intensity components
    """
    # check the format of the spectrum string
    if not re.match(r'^([0-9]+([.][0-9]*)*([eE]-*[0-9]+)*:[0-9]+([.][0-9]*)*([eE]-*[0-9]+)*[ ]*)+$', s):
        msg = f"str_to_ms2: spectrum string not properly formatted: {s}"
        raise ValueError(msg)
    ms_, is_ = [], []
    if ' ' in s:
        for ms2pk in s.split():
            m_, i_ = [float(_) for _ in ms2pk.split(':')]
            ms_.append(m_)
            is_.append(i_)
    else:
        # deal with strings that have only a single mz:intensity pair
        m_, i_ = [float(_) for _ in s.split(':')]
        ms_.append(m_)
        is_.append(i_)
    return np.array([ms_, is_])

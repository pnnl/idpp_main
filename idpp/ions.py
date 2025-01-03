"""
    idpp/ions.py

    Dylan Ross (dylan.ross@pnnl.gov)

    produce ionized structures
"""


import re


def _get_replacements(s, pat, repl):
    """
    yields all possible variations of s with regex pattern (pat)
    replaced by repl (single replacements)
    """
    for mat in re.finditer(pat, s):
        if len(mat.groups()) > 0:
            yield s[:mat.start()] + repl.format(*mat.groups()) + s[mat.end():]
        else:
            yield s[:mat.start()] + repl + s[mat.end():]
        

def _protonate(smi):
    """
    modify a SMILES structure to reflect protonated structure(s)
    use rules to protonate specific groups
    """
    # each rule has a regex pattern and the replacement if it is found
    rules = [
        # acids
        (r'\[O-\]', 'O'),
        # aromatic amines
        (r'([c0-9)])n([c0-9()])', '{}[nH+]{}'),
        # primary amines
        (r'[(]N[)]', '([NH3+])'),
        (r'^N([Cc[S])', '[NH3+]{}'),
        (r'([C0-9])N$', '{}[NH3+]'),
        # secondary amines
        (r'([C)/(])N([Cc[/])', '{}[NH2+]{}'),
        # tertiary amines
        (r'([C)(])N([0-9(])', '{}[NH+]{}'),
        (r'([Cc)])n([0-9(])', '{}[nH+]{}'),
        # if [nH]  or [NH] is ever explicitly included change to H2+
        (r'\[([Nn])H\]', '[{}H2+]'),
    ]
    # convert to set to make sure no duplicate structures
    protomers = set([_ for pat, repl in rules for _ in _get_replacements(smi, pat, repl)])
    return list(protomers) if len(protomers) > 0 else []
    
    
def _ion_adduct(ion):
    """
    modify a SMILES structure to reflect an adduct with another ion
    calling this function with a specific ion returns a function
    that takes a SMILES structure and returns the ionized SMILES structure
    """
    return lambda smi: [smi + '.{}'.format(ion)]
    

def ionize_smi(smi, adduct):
    """
    takes input SMILES structure, attempts computing ionized 
    structures for specified adducts and returns a dictionary 
    mapping adduct to ionized SMILES structure
    
    Parameters
    ----------
    smi : ``str``
        SMILES structure (must be neutral)
    adduct : ``str``
        adduct to compute ionized SMILES structures for
        
    Returns
    -------
    ions : ``list(str)``
        list of ionized SMILES structures
    """
    # individual functions for handling different adducts
    # each should take SMILES structure as an arg and return 
    # either a list of SMILES structure(s) as strings or None
    ionize_funcs = {
        '[M+H]+': _protonate,
        '[M+Na]+': _ion_adduct('[Na+]'),
        '[M+K]+': _ion_adduct('[K+]'),
        '[M+NH4]+': _ion_adduct('[NH4+]'),
        '[M+HCOO]-': _ion_adduct('C(=O)[O-]'),
        '[M+CH3COO]-': _ion_adduct('CC(=O)[O-]')
    }
    # ensure specified adduct have defined functions 
    if adduct not in ionize_funcs:
        return []
        # just return empty list instead of raising an exception
        #msg = 'ionize: adduct {} not defined'.format(adduct)
        #raise ValueError(msg)
    # compute ionized structures
    return ionize_funcs[adduct](smi)
    
"""
    idpp/db/builder/metlin_ccs.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Utilities for adding METLIN-CCS dataset to IdPPdb
    - https://metlin.scripps.edu/docs/METLIN-CCS-03-15-2024.xlsx
"""


import os
import errno
from typing import Tuple

import polars

from idpp.db.util import IdPPdb


def _iter_metlin_ccs(metlin_ccs_file: str, dimer_line_params: Tuple[float, float]):
    """
    iterate over rows from the METLIN-CCS file, yielding data from one row at a time
    """
    fix_adducts = {
        "M+Na]": "[M+Na]+",
        "[M-H]": "[M-H]-",
        "[M+H]": "[M+H]+"
    }
    dimer_m, dimer_b = dimer_line_params
    df = polars.read_excel(metlin_ccs_file)[:, [0, 10, 1, 2, 9, 7]]
    for name, adduct, formula, metlin_id, mz, ccs_avg in df.iter_rows():
        if (fixed_adduct := fix_adducts.get(adduct)) is not None:
            # assign dimers by position relative to dimer line
            # points above are dimers and points below are monomers
            # for now exclude the dimer values
            if ccs_avg <= dimer_m * mz + dimer_b:
                yield name, fixed_adduct, formula, str(metlin_id), mz, 1 if fixed_adduct[-1] == "+" else -1, ccs_avg


def add_metlin_ccs_to_idppdb(db: IdPPdb,
                             metlin_ccs_file: str
                             ) -> None :
    """
    Add downloaded METLIN-CCS (.xlsx) to IdPPdb

    Parameters
    ----------
    db : ``IdPPdb``
        IdPP database interface 
    metlin_ccs_file : ``str``
        path to METLIN-CCS download file (.xlsx)
    """
    print("Adding METLIN-CCS to IdPPdb ...")
    # ensure that the METLIN-CCS download file exists
    if not os.path.isfile(metlin_ccs_file):
        raise FileNotFoundError(errno.ENOENT, 
                                os.strerror(errno.ENOENT), 
                                metlin_ccs_file)
    # add source info
    src_id = db.insert_src("METLIN-CCS", "https://doi.org/10.1038/s41592-023-02078-5")
    # add the CCS data
    for name, adduct, form, metlin_id, mz, z, ccs in _iter_metlin_ccs(metlin_ccs_file, 
                                                                      (0.2692, 121.5385)):
        # add the formula
        form_id = db.insert_form(form) #if form is not None else -1
        # add a compound entry
        cmpd_id = db.insert_cmpd(name, form_id=form_id)
        # add adduct entry
        adduct_id = db.insert_adduct(adduct, cmpd_id, mz, z)
        # add ccs
        ccs_id = db.insert_ccs(ccs, adduct_id, src_id)
        # add external identifier
        db.insert_ext_id(cmpd_id, src_id, metlin_id)
    # add a change log entry
    db.insert_change_log_entry("idpp.db.builder.metlin_ccs.add_metlin_ccs_to_idppdb", 
                               "add METLIN-CCS")
    print("... done")

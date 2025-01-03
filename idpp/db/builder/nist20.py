"""
    idpp/db/builder/nist20.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Utilities for adding MS/MS spectra from NIST20 database
"""


import os
import errno
import sqlite3
import json

import numpy as np

from idpp.db.util import IdPPdb
from idpp.db.builder._util import parse_ce


# NIST20 data selection query
_QRY_SEL = """
SELECT
    Name,
    InChIKey,
    NISTNO,
    Formula,
    PrecursorMZ,
    Precursor_type,
    Collision_energy,
    Instrument_type,
    "M/Z",
    Intensity
FROM 
    spectrum
WHERE
    Spectrum_type IN (
        "MS2", 
        "ms2"
    )
    AND Collision_gas="N2"
    AND Precursor_type IN (
        "[M+H]+",
        "[M-H]-",
        "[M+H-H2O]+",
        "[M+Na]+",
        "[2M+H]+",
        "[2M-H]-",
        "[M+OH]-",
        "[M+H-2H2O]+",
        "[M+H-NH3]+",
        "[M-H-CO2]-",
        "[M+2H]2+",
        "[M-H-H2O]-",
        "[M+NH4]+",
        "[M+Cl]-",
        "[M+CHO2]-",
        "[M+C2H4O2]-",
        "[M-2H]2-",
        "[M+2Na-H]+",
        "[M+2Na]2+",
        "[M+K]+",
        "[3M+H]+",
        "[3M-H]-"
    )
"""


def add_nist20_msms_to_idppdb(db: IdPPdb,
                              nist20_dbf: str
                              ) -> None :
    """
    Add NIST20 MS/MS database (a SQLite database file) to IdPPdb

    Parameters
    ----------
    db : ``IdPPdb``
        IdPP database interface 
    nist20_dbf : ``str``
        path to NIST20 database file
    """
    print("Adding NIST20 MS/MS to IdPPdb ...")
    # ensure that the MoNA-export-Experimental_Spectra.json download file exists
    if not os.path.isfile(nist20_dbf):
        # NOTE: Not sure exactly where this file comes from, and how it is produced from 
        #       whatever sort of file is originally distributed when NIST20 was purchased.
        #       Also, now there is a newer version: NIST23
        raise FileNotFoundError(errno.ENOENT, 
                                os.strerror(errno.ENOENT), 
                                nist20_dbf)
    # connect to NIST20 database
    con = sqlite3.connect(nist20_dbf)
    cur = con.cursor()
    # add source info (singular for mapping external IDs)
    single_src_id = db.insert_src("NIST20", "https://www.sisweb.com/software/ms/nist.htm")
    # add source info (separate by Instrument_type for MS2 spectra references)
    src_ids = {
        _: db.insert_src(f"NIST20_{_[0]}", "https://www.sisweb.com/software/ms/nist.htm", 
                         src_notes=f"Instrument_type='{_}'")
        for _ in ["HCD", "Q-TOF", "IT-FT/ion trap with FTMS"]
    }
    # add spectra to the database
    i = 0
    for name, inchi_key, nistno, formula, mz, adduct, ce_str, inst_type, msms_mz, msms_i in cur.execute(_QRY_SEL):
        # add formula, InChI if provided
        form_id = db.insert_form(formula)  # formula always present
        inchi_id = db.insert_inchi(inchi_key) if inchi_key is not None else -1
        # add compounds entry
        # TODO: For now just take the name from the spectrum table, but in the future will grab
        #       more names from the synon table and associate those as well
        cmpd_id = db.insert_cmpd(name, form_id=form_id, inchi_id=inchi_id)
        # add adduct entry
        z = {'+': 1, '-': -1}.get(adduct[-1], 0)
        add_id = db.insert_adduct(adduct, cmpd_id, mz, z)
        # add MS/MS spectrum
        ce = parse_ce(ce_str)
        msms_mz, msms_i =  np.array([json.loads(msms_mz), json.loads(msms_i)])
        _ = db.insert_ms2(msms_mz, msms_i, add_id, src_ids[inst_type], ms2_ce=ce)
        # add external ID if present
        if nistno is not None:
            db.insert_ext_id(cmpd_id, single_src_id, nistno)
        # only print some info every so often
        i += 1
        if i % 100 == 0:
            print(f"\r\tprocessed {i:6d} entries", end="      ")
    print()
    con.close()
    # add a change log entry
    db.insert_change_log_entry("idpp.db.builder.nist20.add_nist20_msms_to_idppdb", 
                               "add NIST20 MS/MS spectra")
    print("... done")


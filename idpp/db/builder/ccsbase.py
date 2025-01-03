"""
    idpp/db/builder/ccsbase.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Utilities for adding CCSbase experimental CCS database to IdPPdb
    - Repository: https://github.com/dylanhross/c3sdb
    - Follow the instructions to build a local copy of C3S.db (a SQLite3 
        database), this is the input used for adding entries to IdPPdb
"""


import os
import errno
import sqlite3
import json

from idpp.db.util import IdPPdb


_QRY_SEL = """
SELECT 
    g_id,
    name,
    adduct,
    z,
    mz,
    ccs,
    smi,
    src_tag,
    ccs_type,
    ccs_method
FROM 
    master
"""


def _ccsbase_iter(cur):
    """ Query the C3S.db and yield data one row at a time """
    for row in cur.execute(_QRY_SEL):
        yield row


def _make_src(src_tag, ccs_type, ccs_method):
    """ returns src_name, reference, and notes from C3S.db source info """
    src_name = f"CCSbase_{src_tag}"
    src_ref = "https://github.com/dylanhross/c3sdb"
    src_notes = json.dumps({"ccs_type": ccs_type, "ccs_method": ccs_method})
    return src_name, src_ref, src_notes


def add_ccsbase_to_idppdb(db: IdPPdb,
                          c3sdb_file: str
                          ) -> None :
    """
    Add downloaded MoNA experimental MS/MS database (a JSON file) to IdPPdb

    Parameters
    ----------
    db : ``IdPPdb``
        IdPP database interface 
    c3sdb_file : ``str``
        path to local C3S.db database
    """
    print("Adding CCSbase to IdPPdb ...")
    # ensure that the C3S.db file
    if not os.path.isfile(c3sdb_file):
        # NOTE: Not sure exactly where this file comes from, and how it is produced from 
        #       whatever sort of file is originally distributed when NIST20 was purchased.
        #       Also, now there is a newer version: NIST23
        raise FileNotFoundError(errno.ENOENT, 
                                os.strerror(errno.ENOENT), 
                                c3sdb_file)
    # connect to C3S.db database
    con = sqlite3.connect(c3sdb_file)
    cur = con.cursor()
    # add source info (singular for mapping external IDs)
    single_src_id = db.insert_src("CCSbase", "https://github.com/dylanhross/c3sdb")
    # add the CCS data
    for g_id, name, adduct, z, mz, ccs, smi, src_tag, ccs_type, ccs_method in _ccsbase_iter(cur):
        # add source info
        src_name, src_ref, src_notes = _make_src(src_tag, ccs_type, ccs_method)
        src_id = db.insert_src(src_name, src_ref, src_notes=src_notes)
        # add smiles entry
        smi_id = db.insert_smi(smi) if smi is not None else -1
        # add a compound entry
        cmpd_id = db.insert_cmpd(name, smi_id=smi_id)
        # add adduct entry
        adduct_id = db.insert_adduct(adduct, cmpd_id, mz, z)
        # add ccs
        ccs_id = db.insert_ccs(ccs, adduct_id, src_id)
        # add external id
        db.insert_ext_id(cmpd_id, single_src_id, g_id)
    # clean up
    cur.close()
    # add a change log entry
    db.insert_change_log_entry("idpp.db.builder.ccsbase.add_ccsbase_to_idppdb", 
                               "add CCSbase experimental CCS")
    print("... done")


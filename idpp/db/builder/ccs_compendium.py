"""
    idpp/db/builder/ccs_compendium.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Utilities for adding McLean Group CCS Compendium to IdPPdb
    - https://mcleanresearchgroup.shinyapps.io/CCS-Compendium
"""


import csv
import os
import errno

from idpp.db.util import IdPPdb


def _fix_adduct(adduct):
    """ fix some of the adducts that have weird formats """
    if adduct == "M+": 
        return "[M]+"
    else:
        return adduct.replace("(", "[").replace(")", "]")
        

def _compendium_iter(compendium_file: str):
    """ iterate over rows from the compendium file, yielding data from one row at a time """
    with open (compendium_file, "r") as f:
        rdr = csv.reader(f, delimiter=",")
        # consume the header row
        _ = next(rdr)
        # "Compound","Neutral.Formula",_,"InChi","InChiKey",_,"mz",_,"Ion.Species.Agilent","Charge",_,"CCS"
        # ,_,_,_,_,"Kingdom","Super.Class","Class","Subclass",___
        for (name, formula, _, inchi, inchikey, _, mz, _, adduct, z, _, ccs, 
                _, _ , _, _, kingdom, supercls, cls, subcls, *_) in rdr:
            ccs = float(ccs)
            z = int(z)
            mz = float(mz)
            form = formula if formula != "" else None
            adduct = _fix_adduct(adduct)
            cls_info = [f"{k}:{v}" for k, v in zip(["kingdom", "superclass", "class", "subclass"], 
                                                   [kingdom, supercls, cls, subcls]) if v != ""]
            yield name, form, inchi, inchikey, cls_info, mz, adduct, z, ccs
        

def add_ccs_compendium_to_idppdb(db: IdPPdb,
                                 compendium_file: str
                                 ) -> None :
    """
    Add downloaded CCS compendium (a .csv) to IdPPdb

    Parameters
    ----------
    db : ``IdPPdb``
        IdPP database interface 
    compendium_file : ``str``
        path to CCS compendium download file (csv)
    """
    print("Adding CCS Compendium to IdPPdb ...")
    # ensure that the CCS compendium download file exists
    if not os.path.isfile(compendium_file):
        raise FileNotFoundError(errno.ENOENT, 
                                os.strerror(errno.ENOENT), 
                                compendium_file)
    # add source info
    src_id = db.insert_src("UnifiedCCSCompendium", "https://doi.org/10.1039/C8SC04396E")
    # add the CCS data
    for name, form, inchi, inchikey, cls_info, mz, adduct, z, ccs in _compendium_iter(compendium_file):
        # add the formula
        form_id = db.insert_form(form) if form is not None else -1
        # add InChI key
        inchi_id = db.insert_inchi(inchikey, inchi=inchi) if inchikey is not None else -1
        # add a compound entry
        cmpd_id = db.insert_cmpd(name, form_id=form_id, inchi_id=inchi_id)
        # add classification if provided
        for class_label in cls_info:
            cls_id = db.insert_class_definition(class_label)
            db.insert_class_label(cls_id, cmpd_id)
        # add adduct entry
        adduct_id = db.insert_adduct(adduct, cmpd_id, mz, z)
        # add ccs
        ccs_id = db.insert_ccs(ccs, adduct_id, src_id)
    # add a change log entry
    db.insert_change_log_entry("idpp.db.builder.ccs_compendium.add_ccs_compendium_to_idppdb", 
                               "add McLean group CCS compendium")
    print("... done")

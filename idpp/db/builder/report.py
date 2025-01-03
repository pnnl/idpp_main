"""
    idpp/db/builder/report.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Utilities for adding RT values from RepoRT (https://github.com/michaelwitting/RepoRT)
    - clone/download the repository and reference its location for adding its datasets to the IdPPdb
"""


import csv
import json
import os
import errno
from typing import Tuple

from idpp.db.util import IdPPdb


def _create_src_entry(report_dir: str, 
                      dataset: str
                      ) -> Tuple[None] :
    """
    Extract data for an entry in the Sources table from a RepoRT dataset

    Parameters
    ----------
    report_dir : ``str``
        path to the RepoRT repository directory
    dataset : ``str``
        RepoRT dataset

    Returns
    -------
    src_name : ``str``
    src_ref : ``str``
    src_notes : ``str``
        data for an entry in the Sources table
    """
    # make sure the source files exist
    dset_info_file = os.path.join(report_dir, f"raw_data/{dataset}/{dataset}_info.tsv")
    if not os.path.isfile(dset_info_file):
        raise FileNotFoundError(errno.ENOENT, 
                                os.strerror(errno.ENOENT), 
                                dset_info_file)
    dset_metadata_file = os.path.join(report_dir, f"raw_data/{dataset}/{dataset}_metadata.tsv")
    if not os.path.isfile(dset_metadata_file):
        raise FileNotFoundError(errno.ENOENT, 
                                os.strerror(errno.ENOENT), 
                                dset_metadata_file)
    # creating a dictionary to store all of the info we want to keep in src_notes
    notes_info = {"RepoRT": dataset}
    # looking through the rows of X_info.tsv to grab src information
    with open(dset_info_file, "r", encoding="utf8") as f:
        file_src = csv.reader(f, delimiter="\t")
        header = next(file_src)
        header_value = next(file_src)
        # creating a dictionary that includes the source name and reference url/doi
        src_data = {header[i]: header_value[i] for i in range(len(header))}
        src_name = src_data["name"]
        src_ref = src_data["url"]
        # adding the source and method of the dataset to notes_info
        notes_info.update(source = src_data["source"])
        notes_info.update(method = src_data["method.type"])
    # looking through the rows of _metadata.tsv to grab src information
    with open(dset_metadata_file, "r", encoding="utf8") as f1:
        file_meta = csv.reader(f1, delimiter="\t")
        meta_header = next(file_meta)
        meta_value = next(file_meta)
        #creating a dictionary that includes the source name and reference url/doi
        meta_data = {meta_header[i]: meta_value[i] for i in range(len(meta_header))}
        #adding in all of the chromatographic parameters we"re interested in seeing in notes
        notes_info.update(column_name = meta_data["column.name"])
        notes_info.update(column_length = meta_data["column.length"])
        notes_info.update(column_id = meta_data["column.id"])
        notes_info.update(column_particle_size = meta_data["column.particle.size"])
        notes_info.update(column_temperature = meta_data["column.temperature"])
        notes_info.update(column_flowrate = meta_data["column.flowrate"])
        notes_info.update(gradient_start_A = meta_data["gradient.start.A"])
        notes_info.update(gradient_start_B = meta_data["gradient.start.B"])
        notes_info.update(gradient_end_A = meta_data["gradient.end.A"])
        notes_info.update(gradient_end_B = meta_data["gradient.end.B"])
    # create a json object from the dict for src_notes    
    src_notes = json.dumps(notes_info)
    return src_name, src_ref, src_notes


def _choose_smi(smi_iso, smi_can):
    """ 
    choose which SMILES structure to keep between isomeric and cannonical (prefer isomeric) 
    """
    # if neither was provided, then return None
    if smi_iso == "" and smi_can == "":
        return None
    if smi_iso != "":
        if smi_can != "":
            # if both provided, prefer isomeric
            return smi_iso
        else:
            # only isomeric provided, so return that
            return smi_iso
    else:
        # only cannonical must have been provided
        # (since we already checked for both empty)
        # so return cannonical
        return smi_can


def _rtdata_iter(report_dir: str, 
                 dataset: str
                 ):
    """ 
    Generator that yields one row at a time from RepoRT-master/raw_data/XXXX/XXXX_rtdata.tsv

    Parameters
    ----------
    report_dir : ``str``
        path to the RepoRT repository directory
    dataset : ``str``
        RepoRT dataset
    
    Yields
    ------
    name, form, smi, inchi_key, inchi, rt, pubchem_id, hmdb_id, lmaps_id, kegg_id
    TODO: add descriptions ^
    """
    # make sure the rtdata file exist
    dset_rtdata_file = os.path.join(report_dir, f"raw_data/{dataset}/{dataset}_rtdata.tsv")
    if not os.path.isfile(dset_rtdata_file):
        raise FileNotFoundError(errno.ENOENT, 
                                os.strerror(errno.ENOENT), 
                                dset_rtdata_file)
    with open(dset_rtdata_file, "r", encoding="utf8") as f:
        _ = next(f)
        rdr = csv.reader(f, delimiter="\t")
        # columns:
        #   id name formula rt 
        #   pubchem.cid pubchem.smiles.isomeric pubchem.smiles.canonical 
        #   pubchem.inchi pubchem.inchikey 
        #   id.hmdb id.lipidmaps id.kegg comment
        for _, name, form, rt, pubchem_id, pubchem_smi_iso, pubchem_smi_can, \
                inchi, inchi_key, _, hmdb_id, lmaps_id, kegg_id, *_ in rdr:
            rt = float(rt)
            # many of these may not be populated
            # return appropriate defaults (usually None) if that is the case
            smi = _choose_smi(pubchem_smi_iso, pubchem_smi_can)
            inchi = inchi if inchi != "" else None
            inchi_key = inchi_key if inchi_key != "" else None
            pubchem_id = pubchem_id if pubchem_id != "" else None
            hmdb_id = hmdb_id if hmdb_id != "" else None
            lmaps_id = lmaps_id if lmaps_id != "" else None
            kegg_id = kegg_id if kegg_id != "" else None
            yield name, form, smi, inchi_key, inchi, rt, pubchem_id, hmdb_id, lmaps_id, kegg_id


_QRY_SEL_GROUPED_SRCS = """
SELECT 
    src_notes, 
    GROUP_CONCAT(src_id) AS src_ids, 
    GROUP_CONCAT(src_name) AS src_names, 
    COUNT(*) AS cnt 
FROM 
    Sources 
GROUP BY 
    src_notes 
HAVING 
    src_notes LIKE "%RepoRT%"
    AND cnt > 1 
ORDER BY 
    cnt DESC;
"""


_QRY_UPDATE_RTS_SRC_IDS = """
UPDATE
    RTs
SET 
    src_id=:new_src_id
WHERE
    src_id=:old_src_id;
"""


def _combined_src_notes(src_notes, src_ids, src_names, src_count):
    return json.dumps({"src_ids": src_ids, "src_names": src_names, "src_count": src_count, "notes": json.loads(src_notes)})


def _combine_duplicate_report_sources(db: IdPPdb
                                      ) -> None :
    """
    Combine RT datasets from RepoRT that have the same chromatographic method
    into new grouped Sources entries. The original source IDs and descriptions
    are retained in the database but entries in RTs are remapped to the new
    grouped sources.

    Parameters
    ----------
    db : ``IdPPdb``
        IdPP database interface 
    """
    # create new combined sources and reassign source IDs from RTs
    for i, (src_notes, src_ids, src_names, src_count) in enumerate(db.cur.execute(_QRY_SEL_GROUPED_SRCS).fetchall()):
        # create new source notes
        new_src_notes = _combined_src_notes(src_notes, src_ids, src_names, src_count)
        new_src_id = db.insert_src(f"RepoRT_meta_{i + 1}", "DOI=?", src_notes=new_src_notes)
        for old_src_id in [int(_) for _ in src_ids.split(",")]:
            db.cur.execute(_QRY_UPDATE_RTS_SRC_IDS, {"new_src_id": new_src_id, "old_src_id": old_src_id})


def add_report_datasets_to_idppdb(db: IdPPdb,
                                  report_dir: str
                                  ) -> None :
    """
    Add RT datasets from RepoRT (https://github.com/michaelwitting/RepoRT)  to IdPPdb

    Parameters
    ----------
    db : ``IdPPdb``
        IdPP database interface 
    report_dir : ``str``
        path to the RepoRT repository directory
    """
    print("Adding RepoRT datasets to IdPPdb ...")
    # make sure raw_data directory exists
    raw_data_dir = os.path.join(report_dir, "raw_data/")
    if not os.path.isdir(raw_data_dir):
        # file OR DIRECTORY not found
        raise FileNotFoundError(errno.ENOENT, 
                                os.strerror(errno.ENOENT), 
                                raw_data_dir)
    # add in some more sources for mapping external IDs 
    lmaps_src_id = db.insert_src("LipidMAPS", "https://www.lipidmaps.org/databases/lmsd")
    kegg_src_id = db.insert_src("KEGG", "https://www.genome.jp/kegg/")
    # fetch source ids for hmdb and pubchem
    hmdb_src_id = db.insert_src("HMDB", "already present in DB")
    pubchem_src_id = db.insert_src("PubChem", "already present in DB")
    # iterate over datasets (directories within raw_data_dir)
    datasets = sorted([_ for _ in os.listdir(raw_data_dir) if os.path.isdir(os.path.join(raw_data_dir, _))])
    n_datasets = len(datasets)
    for i, dataset in enumerate(datasets):
        # add source info for the dataset
        src_name, src_ref, src_notes = _create_src_entry(report_dir, dataset)
        src_id = db.insert_src(src_name, src_ref, src_notes=src_notes)
        # iterate over rows in a dataset
        for name, form, smi, inchi_key, inchi, rt, pubchem_id, hmdb_id, lmaps_id, kegg_id in _rtdata_iter(report_dir, dataset):
            # add the formula
            form_id = db.insert_form(form)
            # add SMILES structure
            smi_id = db.insert_smi(smi) if smi is not None else -1
            # add InChI key
            inchi_id = db.insert_inchi(inchi_key, inchi=inchi) if inchi_key is not None else -1
            # add a compound entry
            cmpd_id = db.insert_cmpd(name, 
                                    form_id=form_id, smi_id=smi_id, inchi_id=inchi_id)
            # add adduct entry
            # this is a "psuedo" adduct since these RT values do not have
            # an explicit ionization state available
            adduct_id = db.insert_adduct("none", cmpd_id, 0., 0)
            # add the RT
            rt_id = db.insert_rt(rt, adduct_id, src_id)
            # add any external IDs that were provided
            if pubchem_id is not None:
                db.insert_ext_id(cmpd_id, pubchem_src_id, pubchem_id)
            if hmdb_id is not None:
                db.insert_ext_id(cmpd_id, hmdb_src_id, hmdb_id)
            if lmaps_id is not None:
                db.insert_ext_id(cmpd_id, lmaps_src_id, lmaps_id)
            if kegg_id is not None:
                db.insert_ext_id(cmpd_id, kegg_src_id, kegg_id)
        print(f"\r\tadded dataset: {i + 1:4d} / {n_datasets}", end="      ")
    print(f"\n\tcombining sources with the same chromatographic method ... ", end="")
    _combine_duplicate_report_sources(db)
    print("done")
    # add a change log entry
    db.insert_change_log_entry("idpp.db.builder.report.add_report_datasets_to_idppdb", 
                               "add RepoRT datasets")
    print("... done")



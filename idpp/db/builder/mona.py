"""
    idpp/db/builder/mona.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Utilities for adding MS/MS spectra from MoNA (https://mona.fiehnlab.ucdavis.edu/downloads)
    - Download link: https://mona.fiehnlab.ucdavis.edu/rest/downloads/retrieve/19a23fd5-4e06-4122-ae9d-169198ee9794
"""


from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass
import json
import os
import errno
import glob

from idpp.db.util import IdPPdb
from idpp.db.builder._util import parse_ce, str_to_ms2


def mona_chunks_exist(chunk_dir: str
                      ) -> bool :
    """
    Check that chunked MoNA download files exist
    """
    # if the directory exists and there is at least 1
    # of ####.json files in there, then call it good
    return os.path.isdir(chunk_dir) and len(glob.glob("[0-9][0-9][0-9][0-9].json", root_dir=chunk_dir)) > 0


def chunk_mona_json(mona_export_json: str,
                    chunk_dir: str,
                    ) -> None :
    """
    process a download of experimental MS/MS spectra from MoNA into
    chunks of a managable size (~5 MB)

    Parameters
    ----------
    mona_export_json : ``str``
        path to MoNA-export-Experimental_Spectra.json downloaded file
    chunk_dir : ``str``
        directory to save the chunked files into
    """
    # ensure that the .json download file exists
    if not os.path.isfile(mona_export_json):
        # TODO: Offer to download the file from MoNA? It is pretty big (a few GB) but
        #       could be pretty convenient to have that option.
        raise FileNotFoundError(errno.ENOENT, 
                                os.strerror(errno.ENOENT), 
                                mona_export_json)
    # ensure that the chunk directory exists, if not make it
    if not os.path.isdir(chunk_dir):
        os.mkdir(chunk_dir)
    update_freq = 5
    # load the JSON file
    print("loading JSON ...", end=" ")
    entries = None
    with open(mona_export_json, "r", encoding="utf8") as jf:
        entries = json.load(jf)
    print("done")
    # break into chunks with 256 entries (results in ~900 chunks of a few MB in size)
    chunk_entries = 256
    chunk_n = 1
    chunk = []
    for i, entry in enumerate(entries):
            chunk.append(entry)
            if i % chunk_entries == 0:
                with open(os.path.join(chunk_dir, f"{chunk_n:04d}.json"), "w") as jf:
                    json.dump(chunk, jf)
                chunk = []
                chunk_n += 1
                if chunk_n % update_freq == 0:
                    print('\r{:5d} chunks written {}  '.format(chunk_n, '|/--\\'[chunk_n // update_freq % 5]), end='')
    # write the last chunk, whatever is in there
    with open(os.path.join(chunk_dir, f"{chunk_n:04d}.json"), "w") as jf:
        json.dump(chunk, jf)
    print('\r{:5d} chunks written {}  '.format(chunk_n, '|/--\\'[chunk_n // update_freq % 5]), end='')
    print()


@dataclass
class _ParsedMonaRequired:
    """ required data from a parsed MoNA entry """
    names: List[str]
    mz: float
    adduct: str
    spectrum: str


@dataclass
class _ParsedMonaOptional:
    """ optional data from a parsed MoNA entry """
    ce: Optional[int]
    mona_id: Optional[str]
    formula: Optional[str]
    inchi_key: Optional[str]
    inchi: Optional[str]
    smi: Optional[str]
    class_labels: Optional[List[str]]


def _parse_mona_entry(entry: str
                      ) -> Optional[Tuple[_ParsedMonaRequired, _ParsedMonaOptional]] :
    """
    parse an entry from the MoNA experimental spectra database (JSON), return
    a dict with the relevant info from that entry if successful, else return None

    Parameters
    ----------
    entry : ``str``
        an entry from the exported database, a string in JSON format

    Returns
    -------
    parsed : ``(_ParsedMonaRequired, _ParsedMonaOptional)`` or ``None``
        relevant info for the entry as dataclasses (one for required info, the other for 
        optional) or None if parsing was unsuccessful
    """
    metadata = {_["name"]: _.get("value") for _ in entry["metaData"]}
    # only consider MS2 spectra
    if metadata.get("ms level") != "MS2":
        return None
    # get the required info
    pre_mz = None
    if (pmz := metadata.get("precursor m/z")) is not None:
        pre_mz = float(pmz)
    pre_adduct = metadata.get("precursor type")
    try:
        names = [_["name"] for _ in entry["compound"][0]["names"]]
    except:
        # TODO: is there a specific error to catch here? KeyError maybe?
        names = None
    if len(names) == 0:
        names = None
    spectrum = entry.get("spectrum")
    # validate the required info
    if pre_mz is None or pre_adduct is None or names is None or spectrum is None:
        return None
    # get extra info
    cmpd_metadata = {_["name"]: _.get("value") for _ in entry["compound"][0]["names"]}
    mona_id = entry.get("id")
    formula = cmpd_metadata.get("molecular formula")
    inchi_key = entry["compound"][0].get("inchi")
    inchi = entry["compound"][0].get("inchiKey")
    smi = cmpd_metadata.get("SMILES")
    classification = {_["name"]: _["value"] for _ in entry["compound"][0]["classification"]}
    class_info = [f"{k}:{c}" for k in ["kingdom", "superclass", "class", "subclass", "direct parent"] 
                  if (c := classification.get(k)) is not None]
    ce = parse_ce(metadata.get("collision energy"))
    # assemble and return the parsed entry
    parsed = (
        _ParsedMonaRequired(**{
                                "names": names,
                                "mz": pre_mz,
                                "adduct": pre_adduct,
                                "spectrum": spectrum
                            }),
        _ParsedMonaOptional(**{
                                "ce": ce,
                                "mona_id": mona_id,
                                "formula": formula,
                                "inchi_key": inchi_key,
                                "inchi": inchi,
                                "smi": smi,
                                "class_labels": class_info 
                            })
    )
    return parsed


def add_mona_chunks_to_idppdb(db: IdPPdb,
                              chunk_dir: str
                              ) -> None :
    """
    Add downloaded MoNA experimental MS/MS database (a JSON file) to IdPPdb

    Parameters
    ----------
    db : ``IdPPdb``
        IdPP database interface 
    chunk_dir : ``str``
        path to directory with chunks from MoNA-export-Experimental_Spectra.json downloaded file
    """
    print("Adding MoNA experimental MS/MS to IdPPdb ...")
    # ensure that the MoNA-export-Experimental_Spectra.json download file exists
     # make sure chunk_dir is valid and has chunk files in it
    if not mona_chunks_exist(chunk_dir):
        msg = f"add_mona_chunks_to_idppdb: HMDB chunks not found in {chunk_dir}"
        raise RuntimeError(msg)
    # add source info
    src_id = db.insert_src('MoNa', 'https://mona.fiehnlab.ucdavis.edu/')
    # add spectra to the database
    i = 0
    for chunk in sorted(glob.glob("[0-9][0-9][0-9][0-9].json", root_dir=chunk_dir)):
        entries = None
        with open(os.path.join(chunk_dir, chunk), "r") as jf:
            entries = json.load(jf)
        for entry in entries:
            # parse the entry, if it worked then proceed with adding to the database
            if (parsed := _parse_mona_entry(entry)) is not None:
                p_req, p_opt = parsed
                # add formula, SMILES, InChI if provided
                form_id = db.insert_form(p_opt.formula) if p_opt.formula is not None else -1
                smi_id = db.insert_smi(p_opt.smi) if p_opt.smi is not None else -1
                inchi_id = db.insert_inchi(p_opt.inchi_key, inchi=p_opt.inchi) if p_opt.inchi_key is not None else -1
                # add compounds entry
                # TODO: For now just take the first name, but in the future give the whole list to
                #       insert_cmpd and it will handle adding to the proper compound entry and dealing
                #       with cases where there are multiple names
                name = p_req.names[0]
                cmpd_id = db.insert_cmpd(name, form_id=form_id, smi_id=smi_id, inchi_id=inchi_id)
                # add adduct entry
                z = {'+': 1, '-': -1}.get(p_req.adduct[-1], 0)
                add_id = db.insert_adduct(p_req.adduct, cmpd_id, p_req.mz, z)
                # add MS/MS spectrum
                msms_mz, msms_i =  str_to_ms2(p_req.spectrum)
                _ = db.insert_ms2(msms_mz, msms_i, add_id, src_id, ms2_ce=p_opt.ce)
                # add external ID and classification info if present
                if p_opt.mona_id is not None:
                    db.insert_ext_id(cmpd_id, src_id, p_opt.mona_id)
                if p_opt.class_labels is not None:
                    for class_label in p_opt.class_labels:
                        cls_id = db.insert_class_definition(class_label)
                        db.insert_class_label(cls_id, cmpd_id)
            i += 1
        # report progress
        print('\r\t{:8d} entries processed (last chunk: {:8s}) '.format(i, chunk), end='')
    print('\r\t{:8d} entries processed'.format(i), end='')
    print()
    # add a change log entry
    db.insert_change_log_entry("idpp.db.builder.mona.add_mona_chunks_to_idppdb", 
                               "add MoNA experimental MS/MS spectra")
    print("... done")


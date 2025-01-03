"""
    Perform some cleaning steps on the database prior to running 
    probability analysis:

    - [x] Exclude lipids from Compounds using regexes
    - [x] Replace any compound names that are quoted with the unquoted version
        - also remove any leading/trailing spaces from compound names
    - [x] names with semicolons
        - split on ; keep first part
        - may also have "w/o MS2:" in front
            - which may have '[0-9]+, "' in front of that
    - [x] names following some patterns with prefixes
        - come up with a regex matching things like "NCGC012345-01!compound name[IIN-based:...]"
        - match the core part with the compound name
        - the first part might be "MLS" instead of "NCGC" ("MLS" has slightly different pattern)
        - the last part "[IIN-based:...]" is not always there
    - [x] names with slightly different prefixes
        - come up with a regex something like "NCGC[0-9]+-[0-9]+_<formula>_<compound_name>"
        - extract the formula and map that
        - extract the compound name and rename
    - [x] names are just HMDB IDs
        - compounds that are just named with their HMDB ID can be remapped to other entries
          in Compounds table that have the same HMDB ID
    - [x] after making all compound name updates/replacements, re-map cmpd_ids to get rid
      of any duplicate entries that were introduced
    - [x] get rid of duplicate entries in Compounds that have the same SMILES or InChI structures
"""


import argparse
import re
from typing import List

from idpp.db.util import IdPPdb


# general query for selecting all cmpd_ids and cmpd_names
_QRY_SEL_CMPD_ID_AND_NAME = """--sqlite3
SELECT 
    cmpd_id, 
    cmpd_name
FROM 
    Compounds
;"""


_QRY_DELETE_COMPOUNDS = """--sqlite3
DELETE FROM
    Compounds
WHERE
    cmpd_id IN ({})
;"""


def _drop_lipids_by_regex(db: IdPPdb) -> None :
    pat = re.compile(r'(([a-zA-Z]-?)?[0-9]{2,}:[0-9]+(\([A-Za-z0-9,]+\))?[_\/]?)+')
    drop_ids = []
    for cmpd_id, cmpd_name in db.cur.execute(_QRY_SEL_CMPD_ID_AND_NAME).fetchall():
        if pat.search(cmpd_name) is not None:
            drop_ids.append(str(cmpd_id))
    # drop all the drop ids
    db.cur.execute(_QRY_DELETE_COMPOUNDS.format(",".join(drop_ids)))


# select compounds with quoted names
_QRY_SEL_QUOTED = """--sqlite3
SELECT 
    cmpd_id, 
    cmpd_name
FROM 
    Compounds
WHERE 
    cmpd_name LIKE '"%"'
;"""


# query for updating cmpd_name using cmpd_id
_QRY_UPDATE_CMPD_NAME = """--sqlite3
UPDATE 
    Compounds 
SET 
    cmpd_name=:n
WHERE
    cmpd_id=:i
;"""


def _fix_quoted(db: IdPPdb) -> None :
    for cmpd_id, cmpd_name in db.cur.execute(_QRY_SEL_QUOTED).fetchall():
        db.cur.execute(_QRY_UPDATE_CMPD_NAME, 
                       {"n": cmpd_name.strip('"'), "i": cmpd_id})


# select compounds with leading/trailing spaces in name
_QRY_SEL_LEAD_TRAIL = """--sqlite3
SELECT 
    cmpd_id, 
    cmpd_name
FROM 
    Compounds
WHERE 
    cmpd_name LIKE " %"
    OR cmpd_name LIKE "% "
;"""


def _fix_leading_trailing_spaces(db: IdPPdb) -> None :
    for cmpd_id, cmpd_name in db.cur.execute(_QRY_SEL_LEAD_TRAIL).fetchall():
        db.cur.execute(_QRY_UPDATE_CMPD_NAME, 
                       {"n": cmpd_name.strip(), "i": cmpd_id})


# select out the compounds with semicolons in the name
_QRY_SEL_SEMICOLON_NAMES = """--sqlite3
SELECT 
    cmpd_id, 
    cmpd_name
FROM 
    Compounds
WHERE 
     cmpd_name LIKE "%;%"
;"""


def _split_semicolon_names(db: IdPPdb) -> None :
    cmpd_ids_names = db.cur.execute(_QRY_SEL_SEMICOLON_NAMES).fetchall()
    for cmpd_id, cmpd_name in cmpd_ids_names:
        db.cur.execute(_QRY_UPDATE_CMPD_NAME,
                       {"n": cmpd_name.replace("w/o MS2:", "").split(";")[0], 
                        "i": cmpd_id})


def _extract_ncgc_mls_names(db: IdPPdb) -> None :
    # the "NCGC" and "MLS" entries have slightly different patterns
    # search them separately
    ncgc_pat = re.compile(r'NCGC[0-9]+-[0-9]+!(.+)')
    mls_pat = re.compile(r'MLS[0-9]+-[0-9]+!(.+)[0-9]+-[0-9]+-[0-9]+')
    for cmpd_id, cmpd_name in db.cur.execute(_QRY_SEL_CMPD_ID_AND_NAME).fetchall():
        if (mat := ncgc_pat.match(cmpd_name)) is not None:
            if (extracted := re.sub(r' \[IIN-based:.+', "", mat.group(1))) != "":
                db.cur.execute(_QRY_UPDATE_CMPD_NAME, {"n": extracted, "i": cmpd_id})
                #print(f"{cmpd_name} -> {extracted}")
        # check for the other pattern
        elif (mat := mls_pat.match(cmpd_name)) is not None:
            db.cur.execute(_QRY_UPDATE_CMPD_NAME, {"n": mat.group(1), "i": cmpd_id})
            #print(f"{cmpd_name} -> {mat.group(1)}")


# update a form_id for a Compounds entry identified by cmpd_id
_QRY_UPDATE_FORM_ID = """--sqlite3
UPDATE
    Compounds
SET
    form_id=:f
WHERE
    cmpd_id=:i
;"""


def _extract_ncgc_underscore_names(db: IdPPdb) -> None :
    ncgc_us_pat = re.compile(r'NCGC[0-9]+-[0-9]+_([A-Z0-9]+)_(.+)')
    for cmpd_id, cmpd_name in db.cur.execute(_QRY_SEL_CMPD_ID_AND_NAME).fetchall():
        if (mat := ncgc_us_pat.match(cmpd_name)) is not None:
            # update name
            new_name = mat.group(2).split(", ")[0]
            db.cur.execute(_QRY_UPDATE_CMPD_NAME, {"n": new_name, "i": cmpd_id})
            # update formula ID
            form_id = db.insert_form(mat.group(1))
            db.cur.execute(_QRY_UPDATE_FORM_ID, {"f": form_id, "i": cmpd_id})
            #print(f"{cmpd_name} -> {new_name} ({form_id})")


# select pairs of cmpd_ids, the second cmpd_id can be re-mapped to the first cmpd_id
_QRY_SEL_HMDBID_REPLACEMENTS = """--sqlite3
SELECT 
    GROUP_CONCAT(cmpd_id)
FROM 
    ExternalIDs
WHERE
    ext_id LIKE "HMDB%"
GROUP BY
    ext_id
HAVING
    COUNT(*) > 1
;"""


# select a compound name by cmpd_id
_QRY_SEL_CMPD_NAME_BY_ID = """--sqlite3
SELECT 
    cmpd_name
FROM
    Compounds
WHERE
    cmpd_id=?
;""" 


def _associate_hmdb_names(db: IdPPdb) -> None :
    for res1 in db.cur.execute(_QRY_SEL_HMDBID_REPLACEMENTS).fetchall():
        cmpd_id_A, cmpd_id_B = [int(_) for _ in res1[0].split(",")]
        # the compound might have been a lipid that was removed earlier so make sure
        # there was actually an entry to get the name from first
        if len((res2 := db.cur.execute(_QRY_SEL_CMPD_NAME_BY_ID, (cmpd_id_A,)).fetchall())) > 0:
            db.cur.execute(_QRY_UPDATE_CMPD_NAME, {"n": res2[0][0], "i": cmpd_id_B})


# select compound entries with the same names
_QRY_SEL_DUPLICATE_CMPDS = """--sqlite3
SELECT 
    GROUP_CONCAT(cmpd_id), 
    GROUP_CONCAT(form_id), 
    GROUP_CONCAT(smi_id), 
    GROUP_CONCAT(inchi_id) 
FROM 
    Compounds 
GROUP BY 
    cmpd_name 
HAVING 
    COUNT(*) > 1
    AND cmpd_name != ""
;"""


# select compound entries with the same smi_ids
_QRY_SEL_DUPLICATE_CMPDS_BY_SMIS = """--sqlite3
SELECT 
    GROUP_CONCAT(cmpd_id), 
    GROUP_CONCAT(form_id), 
    smi_id, 
    GROUP_CONCAT(inchi_id) 
FROM 
    Compounds 
GROUP BY 
    smi_id
HAVING 
    COUNT(*) > 1
    AND smi_id > 0
;"""


# select compound entries with the same inchi_ids
_QRY_SEL_DUPLICATE_CMPDS_BY_INCHIS = """--sqlite3
SELECT 
    GROUP_CONCAT(cmpd_id), 
    GROUP_CONCAT(form_id), 
    GROUP_CONCAT(smi_id), 
    inchi_id 
FROM 
    Compounds 
GROUP BY 
    inchi_id
HAVING 
    COUNT(*) > 1
    AND inchi_id > 0
;"""


# update formulas, smiles, and inchis identifiers for a specified compound
_QRY_UPDATE_CMPD_IDENTIFIERS  = """--sqlite3
UPDATE
    Compounds 
SET 
    form_id=:f,
    smi_id=:s,
    inchi_id=:i
WHERE 
    cmpd_id=:c
;"""


def _keep_identifier(identifier_group: str) -> int : 
    if len(identifiers := [int(_) for _ in identifier_group.split(",") if int(_) > 0]) > 0:
        return min(identifiers)
    else:
        return -1


# select all Adducts entries that have a specified cmpd_id
_QRY_SEL_ADDUCT_BY_CMPD_ID = """--sqlite3
SELECT
    adduct_id,
    adduct
FROM 
    Adducts
WHERE
    cmpd_id=?
;"""


# query to update the adducts table with a new cmpd_id
_QRY_UPDATE_ADDUCTS = """--sqlite3
UPDATE
    Adducts
SET
    cmpd_id=:cid
WHERE
    adduct_id=:aid
;"""


# queries for updating property tables
_QRY_UPDATE_RTS = """--sqlite3
UPDATE 
    RTs
SET
    adduct_id=:new
WHERE
    adduct_id=:old
;"""

_QRY_UPDATE_CCSS = """--sqlite3
UPDATE 
    CCSs
SET
    adduct_id=:new
WHERE
    adduct_id=:old
;"""

_QRY_UPDATE_MS2S = """--sqlite3
UPDATE 
    MS2Spectra
SET
    adduct_id=:new
WHERE
    adduct_id=:old
;"""


_QRY_DELETE_ADDUCTS = """--sqlite3
DELETE FROM
    Adducts
WHERE
    adduct_id IN ({})
;"""


def _remap_adducts(db: IdPPdb, keep_cmpd_id: int, drop_cmpd_ids: List[int]) -> None :
    # select out the adduct entries to keep 
    # adduct mapped to adduct_id
    keep_adducts = {
        adduct: adduct_id
        for adduct_id, adduct
        in db.cur.execute(_QRY_SEL_ADDUCT_BY_CMPD_ID, (keep_cmpd_id,)).fetchall()
    }
    #print(f"\t(pre) {keep_adducts=}")
    drop_adduct_ids = []
    for drop_cmpd_id in drop_cmpd_ids:
        for drop_adduct_id, drop_adduct in db.cur.execute(_QRY_SEL_ADDUCT_BY_CMPD_ID, 
                                                          (drop_cmpd_id,)).fetchall():
            if (update_adduct_id := keep_adducts.get(drop_adduct)) is not None:
                # propagate down to property tables
                for qry in [_QRY_UPDATE_RTS, _QRY_UPDATE_CCSS, _QRY_UPDATE_MS2S]:
                    db.cur.execute(qry, {"new": update_adduct_id, "old": drop_adduct_id})
                # add to drop list
                drop_adduct_ids.append(drop_adduct_id)
            else:
                # only need to update the compound_id in the adduct entry and also 
                # update keep_adducts dict
                db.cur.execute(_QRY_UPDATE_ADDUCTS, {"cid": keep_cmpd_id, "aid": drop_adduct_id})
                keep_adducts[drop_adduct] = drop_adduct_id
                # do not add to drop list
    # now the adduct_ids in drop list can be safely deleted
    batch_size = 100
    n_drop = len(drop_adduct_ids)
    for i in range(n_drop // batch_size + (1 if n_drop % batch_size else 0)):
        i_start = i * batch_size
        i_end = i_start + batch_size
        i_end = i_end if i_end < n_drop - 1 else i_end + (n_drop % batch_size)
        drop_adduct_ids_batch = drop_adduct_ids[i * batch_size:i_end]
        db.cur.execute(_QRY_DELETE_ADDUCTS.format(
                            ("?," * len(drop_adduct_ids_batch)).rstrip(",")
                    ),
                    drop_adduct_ids_batch)
    #print(f"\t(post) {keep_adducts=}")
    #print(f"\t{drop_adduct_ids=}")


# select all external identifiers using cmpd_id
_QRY_SEL_EXTID = """--sqlite3
SELECT 
    src_id,
    ext_id
FROM
    ExternalIDs
WHERE
    cmpd_id=?
;"""


# delete an ExternalIDs entry using cmpd_id
_QRY_DELETE_EXTID = """--sqlite3
DELETE FROM
    ExternalIDs
WHERE
    cmpd_id=?
;"""


def _remap_extids(db: IdPPdb, keep_cmpd_id: int, drop_cmpd_ids: List[int]) -> None :
    # go through the drop compound ids and just try to add their external ids 
    # into the database using the keep compound id, then it will be safe to drop
    # all of the drop ids from ExternalIDs
    # let the IdPPdb object handle avoiding duplicates
    for drop_cmpd_id in drop_cmpd_ids:
        for src_id, ext_id in db.cur.execute(_QRY_SEL_EXTID, (drop_cmpd_id,)).fetchall():
            # this will add a new entry if one wasnt mapped before otherwise it will do nothing
            db.insert_ext_id(keep_cmpd_id, src_id, ext_id)
        # now safe to delete all the ExternalIDs entries for this cmpd_id
        db.cur.execute(_QRY_DELETE_EXTID, (drop_cmpd_id,))


# delete a ClassLabels entry using cmpd_id
_QRY_DELETE_CLSLBL = """--sqlite3
DELETE FROM
    ClassLabels
WHERE
    cmpd_id=?
;"""


def _remap_clslbls(db: IdPPdb, keep_cmpd_id: int, drop_cmpd_ids: List[int]) -> None :
    # in this case we cannot let the IdPPdb object handle avoiding duplicates
    # and this information is a bit low priority
    # so for now I will just be lazy and directly drop the classification labels
    # for all of the drop_cmpd_ids without trying to save/remap them first
    for drop_cmpd_id in drop_cmpd_ids:
        # remap class labels here...
        # now safe to delete all the ExternalIDs entries for this cmpd_id
        db.cur.execute(_QRY_DELETE_CLSLBL, (drop_cmpd_id,))


def _remap_duplicate_compounds(db: IdPPdb) -> None :
    for cmpd_ids, form_ids, smi_ids, inchi_ids in db.cur.execute(_QRY_SEL_DUPLICATE_CMPDS).fetchall():
        # get the cmpd_id to keep and the rest will be dropped
        keep_cmpd_id, *drop_cmpd_ids = [int(_) for _ in cmpd_ids.split(",")]
        #print(f"{keep_cmpd_id=} {drop_cmpd_ids=}")
        # update the entry for keep_cmpd_id with aggregated identifiers from the group
        db.cur.execute(_QRY_UPDATE_CMPD_IDENTIFIERS, {
            "f": _keep_identifier(form_ids),
            "s": _keep_identifier(smi_ids),
            "i": _keep_identifier(inchi_ids),
            "c": keep_cmpd_id
        })
        # remap Adducts, ExternalIDs, and ClassLabels tables
        _remap_adducts(db, keep_cmpd_id, drop_cmpd_ids)
        _remap_extids(db, keep_cmpd_id, drop_cmpd_ids)
        _remap_clslbls(db, keep_cmpd_id, drop_cmpd_ids)
        # finally, delete the drop compounds
        batch_size = 100
        n_drop = len(drop_cmpd_ids)
        for i in range(n_drop // batch_size + (1 if n_drop % batch_size else 0)):
            i_start = i * batch_size
            i_end = i_start + batch_size
            i_end = i_end if i_end < n_drop - 1 else i_end + (n_drop % batch_size)
            drop_cmpd_ids_batch = drop_cmpd_ids[i * batch_size:i_end]
            db.cur.execute(_QRY_DELETE_COMPOUNDS.format(
                                ("?," * len(drop_cmpd_ids_batch)).rstrip(",")
                        ),
                        drop_cmpd_ids_batch)
            

def _remap_duplicate_compounds_by_smis(db: IdPPdb) -> None :
    for cmpd_ids, form_ids, smi_id, inchi_ids in db.cur.execute(_QRY_SEL_DUPLICATE_CMPDS_BY_SMIS).fetchall():
        # get the cmpd_id to keep and the rest will be dropped
        keep_cmpd_id, *drop_cmpd_ids = [int(_) for _ in cmpd_ids.split(",")]
        #print(f"{keep_cmpd_id=} {drop_cmpd_ids=}")
        # update the entry for keep_cmpd_id with aggregated identifiers from the group
        db.cur.execute(_QRY_UPDATE_CMPD_IDENTIFIERS, {
            "f": _keep_identifier(form_ids),
            "s": smi_id,
            "i": _keep_identifier(inchi_ids),
            "c": keep_cmpd_id
        })
        # remap Adducts, ExternalIDs, and ClassLabels tables
        _remap_adducts(db, keep_cmpd_id, drop_cmpd_ids)
        _remap_extids(db, keep_cmpd_id, drop_cmpd_ids)
        _remap_clslbls(db, keep_cmpd_id, drop_cmpd_ids)
        # finally, delete the drop compounds
        batch_size = 100
        n_drop = len(drop_cmpd_ids)
        for i in range(n_drop // batch_size + (1 if n_drop % batch_size else 0)):
            i_start = i * batch_size
            i_end = i_start + batch_size
            i_end = i_end if i_end < n_drop - 1 else i_end + (n_drop % batch_size)
            drop_cmpd_ids_batch = drop_cmpd_ids[i * batch_size:i_end]
            db.cur.execute(_QRY_DELETE_COMPOUNDS.format(
                                ("?," * len(drop_cmpd_ids_batch)).rstrip(",")
                        ),
                        drop_cmpd_ids_batch)
        

def _remap_duplicate_compounds_by_inchis(db: IdPPdb) -> None :
    for cmpd_ids, form_ids, smi_ids, inchi_id in db.cur.execute(_QRY_SEL_DUPLICATE_CMPDS_BY_INCHIS).fetchall():
        # get the cmpd_id to keep and the rest will be dropped
        keep_cmpd_id, *drop_cmpd_ids = [int(_) for _ in cmpd_ids.split(",")]
        #print(f"{keep_cmpd_id=} {drop_cmpd_ids=}")
        # update the entry for keep_cmpd_id with aggregated identifiers from the group
        db.cur.execute(_QRY_UPDATE_CMPD_IDENTIFIERS, {
            "f": _keep_identifier(form_ids),
            "s": _keep_identifier(smi_ids),
            "i": inchi_id,
            "c": keep_cmpd_id
        })
        # remap Adducts, ExternalIDs, and ClassLabels tables
        _remap_adducts(db, keep_cmpd_id, drop_cmpd_ids)
        _remap_extids(db, keep_cmpd_id, drop_cmpd_ids)
        _remap_clslbls(db, keep_cmpd_id, drop_cmpd_ids)
        # finally, delete the drop compounds
        batch_size = 100
        n_drop = len(drop_cmpd_ids)
        for i in range(n_drop // batch_size + (1 if n_drop % batch_size else 0)):
            i_start = i * batch_size
            i_end = i_start + batch_size
            i_end = i_end if i_end < n_drop - 1 else i_end + (n_drop % batch_size)
            drop_cmpd_ids_batch = drop_cmpd_ids[i * batch_size:i_end]
            db.cur.execute(_QRY_DELETE_COMPOUNDS.format(
                                ("?," * len(drop_cmpd_ids_batch)).rstrip(",")
                        ),
                        drop_cmpd_ids_batch)
            

def _main():
    # get the target database from the command-line
    parser = argparse.ArgumentParser()
    parser.add_argument("idpp_db", help="target IdPP database file")
    parser.add_argument("--dry-run", help="do not commit changes to the database", action="store_true")
    args = parser.parse_args()
    # connect to database
    db = IdPPdb(args.idpp_db, enforce_idpp_ver=False, combine_ms2=False)
    # drop the lipids
    _drop_lipids_by_regex(db)
    print("Dropped lipids.")
    # fix quoted compound names and compound names with leading or trailing spaces
    _fix_quoted(db)
    print("Fixed compounds with quoted names.")
    _fix_leading_trailing_spaces(db)
    print("Fixed compound names with leading or trailing spaces.")
    # split the compound names with semicolons in them
    _split_semicolon_names(db)
    print("Split compound names with semicolons.")
    # extract the names from "NCGC" or "MLS" entries
    _extract_ncgc_mls_names(db)
    print("Extracted NCGC and MLS names.")
    # extract the "NCGC" names with underscores
    _extract_ncgc_underscore_names(db)
    print("Extracted NCGC names with underscores.")
    # fix entries with HMDB ID as a name
    _associate_hmdb_names(db)
    print("Fixed entries with HMDB ID as compound name.")
    # remap duplicate compounds
    _remap_duplicate_compounds(db)
    print("Remapped duplicate compounds.")
    _remap_duplicate_compounds_by_smis(db)
    print("Remapped duplicate compounds (grouped by SMILES structure).")
    _remap_duplicate_compounds_by_inchis(db)
    print("Remapped duplicate compounds (grouped by InChI structure).")
    # add a change log entry (and update version)
    # this script's docstring gives a nice summary of the changes, so just use that for the notes
    db.insert_change_log_entry("clean_db.py", __doc__)
    print("Added change log entry.")
    # clean up
    if not args.dry_run:
        db.commit()
        # reclaim disk space from deleting a bunch of rows
        db.vacuum()
        print("Committed changes.")
    db.close(ignore_uncommitted_changes=True)


if __name__ == "__main__":
    _main()

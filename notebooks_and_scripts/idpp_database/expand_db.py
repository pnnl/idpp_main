"""
"""


import csv
import sqlite3
import re

import numpy as np

from idpp.db.util import IdPPdb


# define queries
_QUERIES = {
    "old": {
        "cmpd_name_from_adduct_id": 
            """--sqlite3
                SELECT
                    cmpd_name,
                    adduct,
                    adduct_mz, 
                    adduct_z
                FROM 
                    Compounds
                    JOIN 
                        Adducts USING(cmpd_id)
                WHERE
                    adduct_id=?
            ;""",
        "cmpd_name_from_adduct_id_2": 
            """--sqlite3
                SELECT
                    cmpd_name
                FROM 
                    Compounds
                    JOIN 
                        Adducts USING(cmpd_id)
                WHERE
                    adduct_id=?
            ;""",
    },
    "new": {
        "cmpd_id_from_cmpd_name": 
            """--sqlite3
                SELECT 
                    cmpd_id
                FROM 
                    Compounds
                WHERE 
                    cmpd_name=?
            ;""",
    }
    # "?": 
    #     """--sqlite3
    #         ?
    #     ;""",
}


# lipid name regex
_LIPID_PAT = re.compile(r'(([a-zA-Z]-?)?[0-9]{2,}:[0-9]+(\([A-Za-z0-9,]+\))?[_\/]?)+')


def _pred_ccs_or_rt_iter(pred_ccs_or_rt_file):
    """ iterator for predicted CCS/RT values """
    with open(pred_ccs_or_rt_file, "r") as f:
        _ = next(f)
        for _ , ccs_or_rt, adduct_id, _ in csv.reader(f):
            if ccs_or_rt != "":
                yield float(ccs_or_rt), int(adduct_id)


def _pred_ms2_iter(pred_ms2_file):
    """ iterator for predicted MS2 spectra """
    with open(pred_ms2_file, "r") as f:
        _ = next(f)
        for _ , msms_mz_str, msms_i_str, adduct_id, _, ce in csv.reader(f):
            yield np.fromstring(msms_mz_str, sep=","), np.fromstring(msms_i_str, sep=","), int(adduct_id), float(ce)


def _name_looks_like_lipid(name):
    return _LIPID_PAT.search(name) is not None


def _add_ccs_predictions(old_db_cur, new_db: IdPPdb, src_name, predicted_file):
    n = 0
    src_id = new_db.insert_src(src_name, "",)
    for ccs, old_adduct_id in _pred_ccs_or_rt_iter(predicted_file):
        name, adduct, mz, z = old_db_cur.execute(_QUERIES["old"]["cmpd_name_from_adduct_id"], 
                                                 (old_adduct_id,)).fetchall()[0]
        # do not bother trying to search compound names that look like lipids
        if not _name_looks_like_lipid(name) and (
            len(cmpd_ids := new_db.cur.execute(_QUERIES["new"]["cmpd_id_from_cmpd_name"],
                                               (name,)).fetchall()) > 0
        ):
            new_cmpd_id = cmpd_ids[0][0]
            # add adduct entry
            adduct_id = new_db.insert_adduct(adduct, new_cmpd_id, mz, int(z))
            # add CCS entry
            ccs_id = new_db.insert_ccs(ccs, adduct_id, src_id)
            # track how many entries were added
            n += 1
        #     print(f"{name=} {old_adduct_id=}", end=" ")
        #     print(f"{new_cmpd_id=}", end=" ") 
        #     print(f"{adduct=}", end=" ")
        #     print()
        # if n >= 100:
        #     break
    return n


def _add_rt_predictions(old_db_cur, new_db: IdPPdb, src_name, predicted_file):
    n = 0
    src_id = new_db.insert_src(src_name, "",)
    for rt, old_adduct_id in _pred_ccs_or_rt_iter(predicted_file):
        name = old_db_cur.execute(_QUERIES["old"]["cmpd_name_from_adduct_id_2"], 
                                  (old_adduct_id,)).fetchall()[0][0]
        # do not bother trying to search compound names that look like lipids
        if not _name_looks_like_lipid(name) and (
            len(cmpd_ids := new_db.cur.execute(_QUERIES["new"]["cmpd_id_from_cmpd_name"],
                                               (name,)).fetchall()) > 0
        ):
            new_cmpd_id = cmpd_ids[0][0]
            # add adduct entry
            adduct_id = new_db.insert_adduct("none", new_cmpd_id, 0., 0)
            # add CCS entry
            rt_id = new_db.insert_rt(rt, adduct_id, src_id)
            # track how many entries were added
            n += 1
        #     print(f"{name=} {old_adduct_id=}", end=" ")
        #     print(f"{new_cmpd_id=}", end=" ") 
        #     print()
        # if n >= 100:
        #     break
    return n


def _add_ms2_predictions(old_db_cur, new_db: IdPPdb, src_name, predicted_file):
    n = 0
    src_id = new_db.insert_src(src_name, "",)
    for msms_mz, msms_i, old_adduct_id, ce in _pred_ms2_iter(predicted_file):
        name, adduct, z, mz = old_db_cur.execute(_QUERIES["old"]["cmpd_name_from_adduct_id"], 
                                                 (old_adduct_id,)).fetchall()[0]
        # do not bother trying to search compound names that look like lipids
        if not _name_looks_like_lipid(name) and (
            len(cmpd_ids := new_db.cur.execute(_QUERIES["new"]["cmpd_id_from_cmpd_name"],
                                               (name,)).fetchall()) > 0
        ):
            new_cmpd_id = cmpd_ids[0][0]
            # add adduct entry
            adduct_id = new_db.insert_adduct(adduct, new_cmpd_id, mz, int(z))
            # add MS2 entry
            ms2_id = new_db.insert_ms2(msms_mz, msms_i, adduct_id, src_id, ms2_ce=ce)
            # track how many entries were added
            n += 1
        #     print(f"{name=} {old_adduct_id=}", end=" ")
        #     print(f"{new_cmpd_id=}", end=" ") 
        #     print(f"{adduct=} {ce=}", end=" ")
        #     print()
        # if n >= 100:
        #     break
    return n


def _main():
    # connect to the old db
    old_db_con = sqlite3.connect("_predicted/idpp.db")
    old_db_cur = old_db_con.cursor()
    # connect to new db
    new_db = IdPPdb("idpp_cleaned_expanded.db")
    # add predictions from CCSbase
    print("CCSbase ... ", end="")
    n = _add_ccs_predictions(old_db_cur, new_db, "PREDICTED_CCSbase", "_predicted/c3sdb.tsv")
    print(f"{n=}")
    # add predictions from DarkChem
    print("DarkChem ... ", end="")
    n = _add_ccs_predictions(old_db_cur, new_db, "PREDICTED_DarkChem", "_predicted/darkchem.tsv")
    print(f"{n=}")
    # add predictions from DeepCCS
    print("DeepCCS ... ", end="")
    n = _add_ccs_predictions(old_db_cur, new_db, "PREDICTED_DeepCCS", "_predicted/deepccs.tsv")
    print(f"{n=}")
    # add predictions from rtp
    print("rtp ... ", end="")
    n = _add_rt_predictions(old_db_cur, new_db, "PREDICTED_idpp_rtp", "_predicted/rtp.tsv")
    print(f"{n=}")
    # add predictions from graff-ms
    print("graff-ms ... ", end="")
    n = _add_ms2_predictions(old_db_cur, new_db, "PREDICTED_graff-ms", "_predicted/graff-ms.tsv")
    print(f"{n=}")
    # clean up
    old_db_con.close()
    new_db.commit()
    new_db.close()


if __name__ == "__main__":
    _main()

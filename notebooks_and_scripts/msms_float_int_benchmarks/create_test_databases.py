"""
"""


import sqlite3
import os
import glob
import json

import numpy as np

from idpp.db.builder.mona import _parse_mona_entry
from idpp.db.builder._util import str_to_ms2


QRY_CREATE_FLOAT = """
CREATE TABLE MS2Fragments (
    ms2_id INT NOT NULL,
    frag_fmz REAL NOT NULL,
    frag_fi REAL NOT NULL
) STRICT;
"""


QRY_CREATE_INT = """
CREATE TABLE MS2Fragments (
    ms2_id INT NOT NULL,
    frag_imz INT NOT NULL,
    frag_ii INT NOT NULL
) STRICT;
"""


# same for both
QRY_INSERT = """ INSERT INTO MS2Fragments VALUES (?,?,?); """


MONA_CHUNK_DIR = (
    "/Users/ross200/Library/CloudStorage/OneDrive-PNNL/"
    "Documents/projects/IdPP/db/integrated_build_system/mona_chunks"
)


def main():
    dbf_float = "bench_float.db"
    dbf_int = "bench_int.db"
    # overwrite if exists
    for dbf in [dbf_float, dbf_int]:
        if os.path.isfile(dbf):
            os.remove(dbf)
    # database connections
    con_float = sqlite3.connect(dbf_float)
    cur_float = con_float.cursor()
    con_int = sqlite3.connect(dbf_int)
    cur_int = con_int.cursor()
    # init databases
    cur_float.execute(QRY_CREATE_FLOAT)
    cur_int.execute(QRY_CREATE_INT)
    # fill databases
    n = 0
    for chunk in sorted(glob.glob("[0-9][0-9][0-9][0-9].json", root_dir=MONA_CHUNK_DIR)):
        with open(os.path.join(MONA_CHUNK_DIR, chunk), "r") as jf:
            for entry in json.load(jf):
                # parse the entry, if it worked then proceed with adding to the database
                if (parsed := _parse_mona_entry(entry)) is not None:
                    n += 1
                    p_req, p_opt = parsed
                    ms2_fmz, ms2_fi = str_to_ms2(p_req.spectrum)
                    ms2_imz = np.round((1e5 * ms2_fmz)).astype(np.int32)
                    ms2_ii = np.round(1e6 * (ms2_fi / sum(ms2_fi))).astype(np.int32)
                    for fmz, fi, imz, ii in zip(ms2_fmz, ms2_fi, ms2_imz, ms2_ii):
                        cur_float.execute(QRY_INSERT, (n, fmz, fi))
                        cur_int.execute(QRY_INSERT, (n, int(imz), int(ii)))
        print('\r\t{:8d} entries processed (last chunk: {:8s}) '.format(n, chunk), end='')
    print('\r\t{:8d} entries processed'.format(n), end='')
    print()
    # clean up
    con_float.commit()
    con_float.close()
    con_int.commit()
    con_int.close()


if __name__ == "__main__":
    main()

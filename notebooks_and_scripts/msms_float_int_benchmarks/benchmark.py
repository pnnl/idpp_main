"""
"""


import glob
import os
import csv
import sqlite3
import time

import numpy as np


# selection queries
QRY_SEL_FLOAT = """ SELECT ms2_id, frag_fmz, frag_fi """
QRY_SEL_INT = """ SELECT ms2_id, frag_imz, frag_ii """


def prep_float_query_spectra():
    spectra = []
    for qtsv in sorted(glob.glob("*.tsv", root_dir="query_spectra")):
        fmz, fi = [], []
        with open(os.path.join("query_spectra", qtsv), "r") as spec_f:
            rdr = csv.reader(spec_f, delimiter="\t")
            for mz, i in rdr:
                fmz.append(float(mz))
                fi.append(float(i))
        spectra.append(np.array([fmz, fi]))
    return spectra


def prep_int_query_spectra():
    spectra = []
    for qtsv in sorted(glob.glob("*.tsv", root_dir="query_spectra")):
        fmz, fi = [], []
        with open(os.path.join("query_spectra", qtsv), "r") as spec_f:
            rdr = csv.reader(spec_f, delimiter="\t")
            for mz, i in rdr:
                fmz.append(int(round(float(mz) * 1e5, 0)))
                fi.append(int(float(i)))
        spectra.append(np.array([fmz, fi]))
    return spectra


def query_database_float(cur_float, query_spectrum, tolerance):
    qfmz, qfi  = query_spectrum
    n = len(qfmz)
    idx = 0
    similarities = {}
    qry = "SELECT ms2_id, frag_fmz, frag_fi FROM MS2Fragments ORDER BY frag_fmz"
    for ms2_id, dbfmz, dbfi in cur_float.execute(qry):
        x = dbfmz - tolerance
        while idx < n and qfmz[idx] < x:
            idx += 1
        if idx == n:
            # stop searching once we have gone all of the way through the query spectrum
            break
        if dbfi == 0.:
            continue
        if qfmz[idx] <= dbfmz - tolerance:
            sum_fi = qfi[idx] + dbfi
            a = sum_fi * np.log2(sum_fi)
            b = qfi[idx] * np.log2(qfi[idx])
            c = dbfi * np.log2(dbfi)
            contribution = a - b - c
            # update the similarities 
            if similarities.get(ms2_id) is not None:
                similarities[ms2_id] += contribution
            else:
                similarities[ms2_id] = contribution
        # else: continue on to next database fragment
    return similarities
            

def query_database_int(cur_float, query_spectrum, tolerance):
    qimz, qii  = query_spectrum
    n = len(qimz)
    idx = 0
    similarities = {}
    qry = "SELECT ms2_id, frag_imz, frag_ii FROM MS2Fragments ORDER BY frag_imz"
    for ms2_id, dbimz, dbii in cur_float.execute(qry):
        dbii = float(dbii)
        x = dbimz - tolerance
        while idx < n and qimz[idx] < x:
            idx += 1
        if idx == n:
            # stop searching once we have gone all of the way through the query spectrum
            break
        if dbii == 0:
            continue
        if qimz[idx] <= dbimz - tolerance:
            sum_fi = qii[idx] + dbii
            a = sum_fi * np.log2(sum_fi)
            b = qii[idx] * np.log2(qii[idx])
            c = dbii * np.log2(dbii)
            contribution = a - b - c
            # update the similarities 
            if similarities.get(ms2_id) is not None:
                similarities[ms2_id] += contribution
            else:
                similarities[ms2_id] = contribution
        # else: continue on to next database fragment
    return similarities


def query_database_float_cache_fragments(cur_float, query_spectra, tolerance):
    t0 = time.time()
    qry = "SELECT ms2_id, frag_fmz, frag_fi FROM MS2Fragments ORDER BY frag_fmz"
    cached = cur_float.execute(qry).fetchall()
    t1 = time.time()
    print(f"pre-load fragments: {t1 - t0:.1f} s")
    for i, (qfmz, qfi)  in enumerate(query_spectra):
        t2 = time.time()
        n = len(qfmz)
        idx = 0
        similarities = {}
        for ms2_id, dbfmz, dbfi in cached:
            x = dbfmz - tolerance
            while idx < n and qfmz[idx] < x:
                idx += 1
            if idx == n:
                # stop searching once we have gone all of the way through the query spectrum
                break
            if dbfi == 0.:
                continue
            if qfmz[idx] <= dbfmz + tolerance:
                sum_fi = qfi[idx] + dbfi
                a = sum_fi * np.log2(sum_fi)
                b = qfi[idx] * np.log2(qfi[idx])
                c = dbfi * np.log2(dbfi)
                contribution = a - b - c
                # update the similarities 
                if similarities.get(ms2_id) is not None:
                    similarities[ms2_id] += contribution
                else:
                    similarities[ms2_id] = contribution
            # else: continue on to next database fragment
        print(f"query spectrum {i + 1}: {time.time() - t2:.1f} s")
        #break
    #print(similarities)
    print(f"total query time: {time.time() - t1:.1f}")


def query_database_int_cache_fragments(cur_float, query_spectra, tolerance):
    t0 = time.time()
    qry = "SELECT ms2_id, frag_imz, frag_ii FROM MS2Fragments ORDER BY frag_imz"
    cached = cur_float.execute(qry).fetchall()
    t1 = time.time()
    print(f"pre-load fragments: {t1 - t0:.1f} s")
    for i, (qimz, qii)  in enumerate(query_spectra):
        t2 = time.time()
        n = len(qimz)
        idx = 0
        similarities = {}
        for ms2_id, dbimz, dbii in cached:
            dbii = float(dbii)
            x = dbimz - tolerance
            while idx < n and qimz[idx] < x:
                idx += 1
            if idx == n:
                # stop searching once we have gone all of the way through the query spectrum
                break
            if dbii == 0:
                continue
            if qimz[idx] <= dbimz + tolerance:
                sum_fi = qii[idx] + dbii
                a = sum_fi * np.log2(sum_fi)
                b = qii[idx] * np.log2(qii[idx])
                c = dbii * np.log2(dbii)
                contribution = a - b - c
                # update the similarities 
                if similarities.get(ms2_id) is not None:
                    similarities[ms2_id] += contribution
                else:
                    similarities[ms2_id] = contribution
            # else: continue on to next database fragment
        print(f"query spectrum {i + 1}: {time.time() - t2:.3f} s")
        #break
    #print(similarities)
    print(f"total query time: {time.time() - t1:.3f}")


def make_int_xlog2x_lookup():
    zero_to_one = np.linspace(0, 2, int(2e6), endpoint=False)
    return np.concatenate([[0], (zero_to_one[1:] * np.log2(zero_to_one[1:]) * 1e6).round(0)]).astype(np.int32)


def query_database_int_cache_fragments_wlookup(cur_float, query_spectra, tolerance):
    t0 = time.time()
    qry = "SELECT ms2_id, frag_imz, frag_ii FROM MS2Fragments ORDER BY frag_imz"
    cached = cur_float.execute(qry).fetchall()
    lookup = make_int_xlog2x_lookup()
    t1 = time.time()
    print(f"pre-load fragments (and make lookup table): {t1 - t0:.1f} s")
    for i, (qimz, qii)  in enumerate(query_spectra):
        t2 = time.time()
        n = len(qimz)
        idx = 0
        similarities = {}
        for ms2_id, dbimz, dbii in cached:
            x = dbimz - tolerance
            while idx < n and qimz[idx] < x:
                idx += 1
            if idx == n:
                # stop searching once we have gone all of the way through the query spectrum
                break
            if dbii == 0:
                continue
            if qimz[idx] <= dbimz + tolerance:
                sum_fi = qii[idx] + dbii
                a = lookup[sum_fi]
                b = lookup[qii[idx]]
                c = lookup[dbii]
                contribution = a - b - c
                # update the similarities 
                if similarities.get(ms2_id) is not None:
                    similarities[ms2_id] += contribution
                else:
                    similarities[ms2_id] = contribution
            # else: continue on to next database fragment
        print(f"query spectrum {i + 1}: {time.time() - t2:.3f} s")
        #break
    #print(similarities)
    print(f"total query time: {time.time() - t1:.3f}")


def main():
    dbf_float = "bench_float.db"
    dbf_int = "bench_int.db"
    # database connections
    con_float = sqlite3.connect(dbf_float)
    cur_float = con_float.cursor()
    con_int = sqlite3.connect(dbf_int)
    cur_int = con_int.cursor()
    # query spectra
    query_spectra_float = prep_float_query_spectra()
    query_spectra_int = prep_int_query_spectra()
    # simulate querying databases and time it
    print("-" * 40)
    print("FLOAT")
    # t0 = time.time()
    # for i, q_spec in enumerate(query_spectra_float):
    #     t1 = time.time()
    #     _ = query_database_float(cur_float, q_spec, 0.05)
    #     print(f"query spectrum {i + 1}: {time.time() - t1:.1f} s")
    # print(f"total time: {time.time() - t0:.1f}")
    query_database_float_cache_fragments(cur_float, query_spectra_float, 0.05)
    print("-" * 40)
    print("INT")
    # t0 = time.time()
    # for i, q_spec in enumerate(query_spectra_int):
    #     t1 = time.time()
    #     _ = query_database_int(cur_int, q_spec, 5000)
    #     print(f"query spectrum {i + 1}: {time.time() - t1:.1f} s")
    # print(f"total time: {time.time() - t0:.1f}")
    query_database_int_cache_fragments(cur_int, query_spectra_int, 5000)
    print("-" * 40)
    print("INT (with lookup)")
    query_database_int_cache_fragments_wlookup(cur_int, query_spectra_int, 5000)
    print("-" * 40)
    # clean up
    con_float.close()
    con_int.close()


if __name__ == "__main__":
    main()

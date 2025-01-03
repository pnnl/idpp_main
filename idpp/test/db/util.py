"""
    idpp/test/db/util.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/db/util.py module
"""


import os
import unittest
from tempfile import TemporaryDirectory, NamedTemporaryFile
import sqlite3

import numpy as np

from idpp import __version__ as IDPP_VER
from idpp.db.util import (
    _MIN_PYTHON_VER,
    _get_tstamp, 
    _db_ver_from_tstamp,
    _add_version_info_and_change_log_entry, 
    create_db, 
    IdPPdb
)


class Test_GetTstamp(unittest.TestCase):
    """ tests for the _get_tstamp function """

    def test_GT_tstamp_correct_format(self):
        """ ensure generated timestamp has correct format """
        # YY/MM/DD-hh:mm
        self.assertRegex(_get_tstamp(), r"[0-9]{2}/[0-9]{2}/[0-9]{2}-[0-9]{2}:[0-9]{2}",
                         msg="timestamp should have format YY/MM/DD-hh:mm")
        

class Test_DbVerFromTstamp(unittest.TestCase):
    """ tests for the _db_ver_from_tstamp function """

    def test_DVFT_correct_format(self):
        """ ensure generated db version has correct format """
        # YYMMDD.HH.MM
        db_ver = _db_ver_from_tstamp(_get_tstamp())
        self.assertRegex(db_ver, r"[0-9]{6}[.][0-9]{2}[.][0-9]{2}",
                         msg=f"db version ({db_ver}) had wrong format")
        
    def test_DVFT_correct_value(self):
        """ ensure generated db version has correct value """
        # YYMMDD.HH.MM
        db_ver = _db_ver_from_tstamp("24/02/28-12:34")
        self.assertEqual(db_ver, "240228.12.34")


class Test_AddVersionInfoAndChangeLogEntry(unittest.TestCase):
    """ tests for the _add_version_info_and_change_log_entry function """

    def test_AVIACLE_correct_version_info_and_change_log(self):
        """ make sure the correct info gets added to an in-memory database """
        # set up a mock database with the expected tables (VersionInfo and ChangeLog)
        db = sqlite3.connect(":memory:")
        cur = db.cursor()
        cur.execute(("CREATE TABLE VersionInfo "
                     "(python_ver TEXT NOT NULL, "
                     "idpp_ver TEXT NOT NULL, "
                     "db_ver TEXT NOT NULL);"))
        cur.execute(("CREATE TABLE ChangeLog "
                     "(tstamp TEXT NOT NULL, "
                     "author TEXT NOT NULL, "
                     "notes TEXT NOT NULL);"))
        # add the version info and change log entry
        _add_version_info_and_change_log_entry(cur)
        # ensure the expected values got added to VersionInfo
        py_ver, idpp_ver, db_ver = cur.execute("SELECT python_ver, idpp_ver, db_ver FROM VersionInfo").fetchall()[0]
        self.assertEqual(py_ver, _MIN_PYTHON_VER, 
                         msg=f"expected Python version: {_MIN_PYTHON_VER}, got Python version: {py_ver}")
        self.assertEqual(idpp_ver, IDPP_VER, 
                         msg=f"expected idpp version: {IDPP_VER}, got idpp version: {idpp_ver}")
        self.assertRegex(db_ver, r"[0-9]{6}[.][0-9]{2}[.][0-9]{2}",
                         msg=f"db version ({db_ver}) had wrong format")
        # ensure the expected values got added to ChangeLog
        tstamp, author, notes = cur.execute("SELECT tstamp, author, notes FROM ChangeLog").fetchall()[0]
        self.assertRegex(tstamp, r"[0-9]{2}/[0-9]{2}/[0-9]{2}-[0-9]{2}:[0-9]{2}",
                         msg="timestamp should have format YY/MM/DD-hh:mm")
        self.assertEqual(author, "idpp.db.util.create_db", 
                         msg="expected author='idpp.db.util.create_db'")
        self.assertEqual(notes, "create database", 
                         msg="expected notes='create database'")
        # clean up
        db.close()
        

class TestCreateDb(unittest.TestCase):
    """ tests for the create_db function """

    def test_CD_overwrite(self):
        """ make sure overwrite flag works as expected """
        # when overwrite flag is False (the default) creating a database using an
        # existing file should raise an error
        with self.assertRaises(RuntimeError,
                               msg="trying to create a database when file exists should raise an error"):
            with NamedTemporaryFile(delete_on_close=False) as ntf:
                # raises an error
                create_db(ntf.name)
        # when overwrite flag is True there is no error
        with NamedTemporaryFile(delete_on_close=False) as ntf:
            # no error
            create_db(ntf.name, overwrite=True)

    # expected tables in the database
    _exp_tables = [
        "_TableDescriptions",
        "VersionInfo",
        "ChangeLog",
        "Sources",
        "Compounds",
        "ClassDefs", 
        "ClassLabels",
        "Adducts",
        "Formulas",
        "ExternalIDs",
        "Smiles",
        "InChIs",
        "RTs",
        "CCSs",
        "MS2Spectra",
        "MS2Fragments",
        "MS2Sources",
        "Datasets",
        "AnalysisResults"
    ]

    def test_CD_expected_db_tables(self):
        """ make sure all of the expected database tables are present """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            dbf = os.path.join(tmp_dir, "test.db")
            # create the database
            create_db(dbf)
            # connect to database
            con = sqlite3.connect(dbf)
            cur = con.cursor()
            qry = "SELECT name FROM sqlite_master WHERE type='table';"
            tables = [_[0] for _ in cur.execute(qry).fetchall()]
            # make sure all expected tables are present
            for exp_table in self._exp_tables:
                self.assertIn(exp_table, tables,
                              msg=f"expected table: {exp_table} not in database")
            # make sure no unexpected tables are present
            for table in tables:
                self.assertIn(table, self._exp_tables,
                              msg=f"unexpected table: {table} in database")
            # clean up
            con.close()

    def test_CD_expected_db_tables_overwrite(self):
        """ make sure all of the expected database tables are present when overwriting """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            dbf = os.path.join(tmp_dir, "test.db")
            # create the database
            create_db(dbf)
            # overwrite the database
            create_db(dbf, overwrite=True)
            # connect to database
            con = sqlite3.connect(dbf)
            cur = con.cursor()
            qry = "SELECT name FROM sqlite_master WHERE type='table';"
            tables = [_[0] for _ in cur.execute(qry).fetchall()]
            # make sure all expected tables are present
            for exp_table in self._exp_tables:
                self.assertIn(exp_table, tables,
                              msg=f"expected table: {exp_table} not in database")
            # make sure no unexpected tables are present
            for table in tables:
                self.assertIn(table, self._exp_tables,
                              msg=f"unexpected table: {table} in database")
            # clean up
            con.close()

    def _check_ver_info_and_change_log(self, cur):
        # ensure the expected values got added to VersionInfo
        py_ver, idpp_ver, db_ver = cur.execute("SELECT python_ver, idpp_ver, db_ver FROM VersionInfo").fetchall()[0]
        self.assertEqual(py_ver, _MIN_PYTHON_VER, 
                         msg=f"expected Python version: {_MIN_PYTHON_VER}, got Python version: {py_ver}")
        self.assertEqual(idpp_ver, IDPP_VER, 
                         msg=f"expected idpp version: {IDPP_VER}, got idpp version: {idpp_ver}")
        self.assertRegex(db_ver, r"[0-9]{6}[.][0-9]{2}[.][0-9]{2}",
                         msg=f"db version ({db_ver}) had wrong format")
        # ensure the expected values got added to ChangeLog
        tstamp, author, notes = cur.execute("SELECT tstamp, author, notes FROM ChangeLog").fetchall()[0]
        self.assertRegex(tstamp, r"[0-9]{2}/[0-9]{2}/[0-9]{2}-[0-9]{2}:[0-9]{2}",
                         msg="timestamp should have format YY/MM/DD-hh:mm")
        self.assertEqual(author, "idpp.db.util.create_db", 
                         msg="expected author='idpp.db.util.create_db'")
        self.assertEqual(notes, "create database", 
                         msg="expected notes='create database'")

    def test_CD_correct_version_info_and_change_log(self):
        """ create a new database and make sure the proper version info and changelog are created """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            dbf = os.path.join(tmp_dir, "test.db")
            # create the database
            create_db(dbf)
            # connect to database
            con = sqlite3.connect(dbf)
            cur = con.cursor()
            self._check_ver_info_and_change_log(cur)
            con.close()

    def test_CD_correct_version_info_and_change_log_overwrite(self):
        """ make sure the proper version info and changelog are created when overwriting an existing file """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            dbf = os.path.join(tmp_dir, "test.db")
            # create the database
            create_db(dbf)
            # overwrite the database
            create_db(dbf, overwrite=True)
            # connect to database
            con = sqlite3.connect(dbf)
            cur = con.cursor()
            self._check_ver_info_and_change_log(cur)
            con.close()


class TestIdPPdb_Smiles(unittest.TestCase):
    """ tests for the IdPPdb class, related to dealing with SMILES structures """

    def test_IDPPDB_Smiles_insert_smi(self):
        """ test inserting SMILES structures into the database -> insert_smi method """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # init database
            dbf = os.path.join(tmp_dir, "test.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            # insert a couple of SMILES structures, they should return the same smi_id
            smis = [
                # tetra-alkyl quaternary ammonium
                "C[N+](CCC)(CCC)CCC",
                "CCC[N+](C)(CCC)CCC",
                "CCC[N+](CCC)(C)CCC",
                "CCC[N+](CCC)(CCC)C", 
            ]
            # insert one entry so that the expected smi_id for the rest is 1
            _ = db.insert_smi("C")
            for smi in smis:
                # all of these should have the same smi_id=1
                self.assertEqual(db.insert_smi(smi), 1)
            # try inserting a malformed SMILES structure 
            # which will return the placeholder id -1
            self.assertEqual(db.insert_smi("not so good"), -1)


class TestIdPPdb_Inchis(unittest.TestCase):
    """ tests for the IdPPdb class, related to dealing with InChI keys/structures """

    def test_IDPPDB_Smiles_insert_smi(self):
        """ test inserting InChI structures into the database -> insert_inchi method """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # init database
            dbf = os.path.join(tmp_dir, "test.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            # if only the inchi_key is provided, the entry gets added as usual
            self.assertEqual(db.insert_inchi("BRMWTNUJHUMWMS-LURJTMIESA-N"), 0)
            # no adding duplicates, should get the same inchi_id back
            self.assertEqual(db.insert_inchi("BRMWTNUJHUMWMS-LURJTMIESA-N"), 0)
            # if given inchi as well but it does not work for whatever reason, the
            # placeholder inchi_id (-1) gets returned
            self.assertEqual(db.insert_inchi("XFNJVJPLKCPIBV-UHFFFAOYSA-N", inchi="bad structure"), -1)
            # inchi_key and valid inchi are provided, should get added as usual
            self.assertEqual(db.insert_inchi("XFNJVJPLKCPIBV-UHFFFAOYSA-N",
                                             "InChI=1S/C3H10N2/c4-2-1-3-5/h1-5H2"), 1)
            # valid inchi is provided but the inchi_key is bad,
            # regenerate the inchi_key then add as usual
            # in this case the id is the same as the previous to avoid duplicates
            self.assertEqual(db.insert_inchi("bad inchi key",
                                             "InChI=1S/C3H10N2/c4-2-1-3-5/h1-5H2"), 1)


class TestIdPPdb_Extids(unittest.TestCase):
    """ tests for the IdPPdb class, related to dealing with external IDs """

    def test_IDPPDB_Extids_insert_extid(self):
        """ test inserting external identifiers into the database -> insert_ext_id method """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # init database
            dbf = os.path.join(tmp_dir, "test.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            # insert a few entries
            entries = [
                (1, 0, "A"),
                (2, 0, "B"),
                (3, 0, "C"),
                (4, 0, "D"),
                (1, 1, "E"),
                (1, 1, "F"),
                (3, 0, "G"),
                (4, 0, "D"),  # this is a duplicate, it should not be added
            ]
            for cmpd_id, src_id, ext_id in entries:
                db.insert_ext_id(cmpd_id, src_id, ext_id)
            # make sure the correct number of entries were added to the database
            self.assertEqual(len(db.cur.execute("SELECT * FROM ExternalIDs").fetchall()), 7)


class TestIdPPdb_MSMS(unittest.TestCase):
    """ tests for the IdPPdb class, related to dealing with MS/MS spectra """

    def _array_elements_equal(self, a, a_exp):
        """ make sure the elements of array a match array a_exp """
        for i, (v, v_exp) in enumerate(zip(a, a_exp)):
            self.assertAlmostEqual(v, v_exp, places=5,
                                   msg=f"element {i} of array ({v}) does not match expected value ({v_exp})")

    def _gen_spectrum(self):
        """ helper method for generating a mock MS/MS spectrum """
        # it can have anywhere from 1 to 100 peaks
        n_peaks = np.random.randint(1, 101)
        mz = np.random.choice(np.arange(50, 1000, 0.05), n_peaks)
        i = 100. * np.random.random(n_peaks)
        idx = np.argsort(mz)
        return mz[idx], i[idx]

    def test_IDPPDB_MSMS_convert_spectrum_to_int_format(self):
        """ test the _convert_spectrum_to_int_format method """
        expected = [
            (
                np.array([[123.45678, 234.56789, 345.67890], [10, 20, 70]]), 
                (
                    [12345678, 23456789, 34567890], 
                    [100000, 200000, 700000]
                )
            ),
            (
                np.array([[12.345678, 23.456789, 34.567890], [10, 20, 10]]), 
                (
                    [1234568, 2345679, 3456789], 
                    [250000, 500000, 250000]
                )
            ),
        ]
        with NamedTemporaryFile(delete_on_close=False) as ntf:
            dbf = ntf.name
            create_db(dbf, overwrite=True)
            db = IdPPdb(dbf)
            for (msms_mz, msms_i), (exp_msms_imz, exp_msms_ii) in expected:
                msms_imz, msms_ii = db._convert_spectrum_to_int_format(msms_mz, msms_i)
                self._array_elements_equal(msms_imz, exp_msms_imz)
                self._array_elements_equal(msms_ii, exp_msms_ii)

    def test_IDPPDB_MSMS_spectrum_int_format_intensity_sum(self):
        """ test that the intensities for integer format MS2 spectra always sum to about 1e6 """
        with NamedTemporaryFile(delete_on_close=False) as ntf:
            dbf = ntf.name
            create_db(dbf, overwrite=True)
            db = IdPPdb(dbf)
            for _ in range(100):
                msms_mz, msms_i = self._gen_spectrum()
                _, msms_ii = db._convert_spectrum_to_int_format(msms_mz, msms_i)
                sum_ii = np.sum(msms_ii)
                # can be al little more than 1000000 due to round-off
                self.assertLessEqual(sum_ii, 1000010)
                # can be a little less than 1000000 due to round-off
                self.assertGreaterEqual(sum_ii, 999990)

    def test_IDPPDB_MSMS_spectrum_int_format_intensity_no_zeros(self):
        """ test that the intensities for integer format MS2 spectra never have zeros """
        with NamedTemporaryFile(delete_on_close=False) as ntf:
            dbf = ntf.name
            create_db(dbf, overwrite=True)
            db = IdPPdb(dbf)
            for _ in range(100):
                msms_mz, msms_i = self._gen_spectrum()
                _, msms_ii = db._convert_spectrum_to_int_format(msms_mz, msms_i)
                for ii in msms_ii:
                    self.assertGreater(ii, 0, msg=f"{msms_i} {msms_ii}")

    def test_IDPPDB_MSMS_convert_spectrum_from_int_format(self):
        """ test the _convert_spectrum_from_int_format method """
        expected = [
            (
                (
                    [12345678, 23456789, 34567890], 
                    [100000, 200000, 700000]
                ),
                np.array([[123.45678, 234.56789, 345.67890], [0.1, 0.2, 0.7]])
            ),
            (
                (
                    [1234568, 2345679, 3456789], 
                    [250000, 500000, 250000]
                ),
                np.array([[12.345678, 23.456789, 34.567890], [0.25, 0.5, 0.25]])
            ),
        ]
        with NamedTemporaryFile(delete_on_close=False) as ntf:
            dbf = ntf.name
            create_db(dbf, overwrite=True)
            db = IdPPdb(dbf)
            for (msms_imz, msms_ii), (exp_msms_mz, exp_msms_i) in expected:
                msms_imz, msms_ii = db._convert_spectrum_from_int_format(msms_imz, msms_ii)
                self._array_elements_equal(msms_imz, exp_msms_mz)
                self._array_elements_equal(msms_ii, exp_msms_i)

    def test_IDPPDB_MSMS_insert_ms2(self):
        """ test inserting MS2 spectra into the database -> insert_ms2 method """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # init database
            dbf = os.path.join(tmp_dir, "test.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            # insert a bunch of spectra, there should at least be no errors
            # print()
            for i in range(100):
                # print("* " * 20)
                msms_mz, msms_i = self._gen_spectrum()
                # print("msms_mz", msms_mz)
                # print("msms_i:", msms_i)
                db.insert_ms2(msms_mz, msms_i, i + 1, -1)
            # TODO: validate some database contents after adding the spectra

    def test_IDPPDB_MSMS_fetch_ms2_data(self):
        """ test fetching MS2 spectra from the database -> fetch_ms2_data method """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # init database
            dbf = os.path.join(tmp_dir, "test.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            # add some compounds/adducts/sources and mass spectra to the DB
            src_id = db.insert_src("source", "none")
            for name in ["compound_A", "compound_B"]:
                cmpd_id = db.insert_cmpd(name)
                for adduct, z in [("[M+H]+", 1), ("[M-H]-", -1)]:
                    adduct_id = db.insert_adduct(adduct, cmpd_id, 420.696, z)
                    for i in range(10):
                        msms_mz, msms_i = self._gen_spectrum()
                        #print(msms_mz, msms_i)
                        ce = 10 * i if (i % 3 == 0 and i > 0) else None
                        src_id = i // 2 if i > 2 else 0
                        db.insert_ms2(msms_mz, msms_i, adduct_id, src_id, ms2_ce=ce)
            db.commit()
            #shutil.copy(dbf, "/Users/ross200/Desktop/test.db")
            # TODO: I have used the above ^ to manually inspect the temporary DB file
            #       and it looks fine to me, still worth doing some more explicit 
            #       validation on the values you get back out from fetching as mentioned
            #       below
            # fetch the MSMS data, there should at least be no errors
            # print("\n", "* " * 20)
            for _ in db.fetch_ms2_data(3):
                # print(db.last_qry)
                # print("--" * 10)
                # print(_)
                # TODO: Validate the batches that are yielded: Shape? Values?
                pass


class TestIdPPdb_ProbabilityAnalysis(unittest.TestCase):
    """ tests for the IdPPdb class, related to probability analysis datasets and results """

    def test_IDPPDB_PA_insert_dataset(self):
        """ test inserting a Dataset entry into the database -> insert_dataset method """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # init database
            dbf = os.path.join(tmp_dir, "test.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            # insert a couple datasets, there should at least be no errors
            dsid_a = db.insert_dataset("a desription of a dataset", 
                                       "SELECT something FROM SomeTable", 
                                       69)
            dsid_b = db.insert_dataset("a desription of a different dataset", 
                                       "SELECT something_else FROM SomeOtherTable", 
                                       420)
            dsid_c = db.insert_dataset("a desription of yet another dataset", 
                                       "SELECT another FROM AnotherTable", 
                                       69420)
            # each dataset should have gotten a unique dataset_id
            self.assertNotEqual(dsid_a, dsid_b)
            self.assertNotEqual(dsid_b, dsid_c)
            self.assertNotEqual(dsid_a, dsid_c)
            # TODO: validate some database contents after adding the datasets

    def test_IDPPDB_PA_insert_analysis_results(self):
        """ test inserting an AnalysisResults entry into the database -> insert_analysis_result method """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # init database
            dbf = os.path.join(tmp_dir, "test.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            # insert a couple results, there should at least be no errors
            counts = np.array([8, 354, 5, 35, 54, 54, 51, 5, 2], dtype=np.int32)
            db.insert_analysis_result(69, counts, 5.)
            db.insert_analysis_result(69, counts, 5., ccs_tol=1.)
            db.insert_analysis_result(69, counts, 5., rt_tol=0.5)
            db.insert_analysis_result(69, counts, 5., ms2_tol=0.8)
            db.insert_analysis_result(69, counts, 5., ccs_tol=1., rt_tol=0.5)
            db.insert_analysis_result(69, counts, 5., ccs_tol=1., ms2_tol=0.8)
            db.insert_analysis_result(69, counts, 5., rt_tol=0.5, ms2_tol=0.8)
            db.insert_analysis_result(69, counts, 5., ccs_tol=1., rt_tol=0.5, ms2_tol=0.8)
            # TODO: validate some database contents after adding the analysis results


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsUtil = unittest.TestSuite()
AllTestsUtil.addTests([
    _loader.loadTestsFromTestCase(Test_GetTstamp),
    _loader.loadTestsFromTestCase(Test_DbVerFromTstamp),
    _loader.loadTestsFromTestCase(Test_AddVersionInfoAndChangeLogEntry),
    _loader.loadTestsFromTestCase(TestCreateDb),
    _loader.loadTestsFromTestCase(TestIdPPdb_Smiles),
    _loader.loadTestsFromTestCase(TestIdPPdb_Inchis),
    _loader.loadTestsFromTestCase(TestIdPPdb_Extids),
    _loader.loadTestsFromTestCase(TestIdPPdb_MSMS),
    _loader.loadTestsFromTestCase(TestIdPPdb_ProbabilityAnalysis),
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsUtil)

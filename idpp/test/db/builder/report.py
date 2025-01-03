"""
    idpp/test/db/builder/report.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/db/builder/report.py module
"""


import unittest
import os
import pathlib
from tempfile import TemporaryDirectory
import contextlib
import io

from idpp.test.__include import TEST_INCLUDE_DIR
from idpp.db.util import create_db, IdPPdb
from idpp.db.builder.report import (
    _create_src_entry, _choose_smi, _rtdata_iter, 
    add_report_datasets_to_idppdb
)


# define path to mock RepoRT repository
_MOCK_REPORT_DIR = os.path.join(TEST_INCLUDE_DIR, "mock_RepoRT-master/")


class Test_MockRepoRtRepo(unittest.TestCase):
    """ tests for built-in mock_RepoRT-master repository """

    def test_MRRR_dir_exists(self):
        """ built-in mock_RepoRT-master directory exists """
        self.assertTrue(os.path.isdir(_MOCK_REPORT_DIR))
    
    # TODO:
    # def test_MRRR_contents_exist(self):
    #     """ check that the expected contents of mock_RepoRT-master exist """
        

class Test_CreateSrcEntry(unittest.TestCase):
    """ tests for the _create_src_entry function """

    def test_CSE_info_or_metadata_files_not_found(self):
        """ should raise an error if the dataset info or metadata files not found """
        # dataset info file exists, metadata does not
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                self.assertRaises(FileNotFoundError):
            os.makedirs(os.path.join(tmp_dir, "raw_data/test/"))
            pathlib.Path.touch(os.path.join(tmp_dir, "raw_data/test/test_info.tsv"))
            _ = _create_src_entry(tmp_dir, "test")
        # dataset metadata file exists, info does not
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                self.assertRaises(FileNotFoundError):
            os.makedirs(os.path.join(tmp_dir, "raw_data/test/"))
            pathlib.Path.touch(os.path.join(tmp_dir, "raw_data/test/test_metadata.tsv"))
            _ = _create_src_entry(tmp_dir, "test")
        # Try the following to test that the temporary files in the above examples get
        # created, it should cause some other error not the FileNotFoundError
        # with TemporaryDirectory() as tmp_dir:
        #     os.makedirs(os.path.join(tmp_dir, "raw_data/test/"))
        #     pathlib.Path.touch(os.path.join(tmp_dir, "raw_data/test/test_info.tsv"))
        #     pathlib.Path.touch(os.path.join(tmp_dir, "raw_data/test/test_metadata.tsv"))
        #     _ = _create_src_entry(tmp_dir, "test")

    def test_CSE_mock_data_noerrs(self):
        """ create source entries for all of the mock data, no errors """
        for i in range(1, 9):
            dataset = f"{int(2**i):04d}"
            src_name, src_ref, src_notes = _create_src_entry(_MOCK_REPORT_DIR, dataset)
            # TODO: furter validate the information returned?


class Test_ChooseSmi(unittest.TestCase):
    """ tests for the _choose_smi function """

    def test_CS_both_empty(self):
        """ isomeric and canonnical SMILES empty """
        self.assertIsNone(_choose_smi("", ""))

    def test_CS_isomeric_empty(self):
        """ isomeric SMILES empty """
        self.assertEqual(_choose_smi("", "ICPOOP"), "ICPOOP")

    def test_CS_isomeric_empty(self):
        """ canonnical SMILES empty """
        self.assertEqual(_choose_smi("POOPIC", ""), "POOPIC")
    
    def test_CS_neither_empty(self):
        """ neither isomeric nor canonnical SMILES empty """
        self.assertEqual(_choose_smi("POOPIC", "ICPOOP"), "POOPIC")


class Test_RtdataIter(unittest.TestCase):
    """ tests for the _rtdata_iter function """

    def test_RI_rtdata_file_not_found(self):
        """ should raise an error if the dataset rtdata file isnt found """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, self.assertRaises(FileNotFoundError):
            os.makedirs(os.path.join(tmp_dir, "raw_data/test/"))
            for _ in _rtdata_iter(tmp_dir, "test"):
                pass
    
    def test_RI_mock_data_noerrs(self):
        """ create iterators for all of the mock datasets, no errors """
        for i in range(1, 9):
            dataset = f"{int(2**i):04d}"
            for _ in _rtdata_iter(_MOCK_REPORT_DIR, dataset):
                pass
            # TODO: furter validate the information returned?


class TestAddReportDatasetsToIdppdb(unittest.TestCase):
    """ tests for the add_report_datasets_to_idppdb function """

    def test_ARDTI_rtdata_file_not_found(self):
        """ should raise an error if the raw_data directory is not found """
        with self.assertRaises(FileNotFoundError):
            with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                    contextlib.redirect_stdout(io.StringIO()) as _:
                dbf = os.path.join(tmp_dir, "idpp.db")
                create_db(dbf)
                db = IdPPdb(dbf)
                add_report_datasets_to_idppdb(db, "this dir does not exist")
    
    def test_ARDTI_add_mock_data(self):
        """ add mock RepoRT datasets to IdPPdb """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                contextlib.redirect_stdout(io.StringIO()) as _:
            # init the database
            dbf = os.path.join(tmp_dir, "idpp.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            add_report_datasets_to_idppdb(db, _MOCK_REPORT_DIR)
            # TODO: Check the counts of some of the tables?


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsReport = unittest.TestSuite()
AllTestsReport.addTests([
    _loader.loadTestsFromTestCase(Test_MockRepoRtRepo),
    _loader.loadTestsFromTestCase(Test_CreateSrcEntry),
    _loader.loadTestsFromTestCase(Test_ChooseSmi),
    _loader.loadTestsFromTestCase(Test_RtdataIter),
    _loader.loadTestsFromTestCase(TestAddReportDatasetsToIdppdb),
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsReport)

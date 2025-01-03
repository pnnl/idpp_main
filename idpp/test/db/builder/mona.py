"""
    idpp/test/db/builder/mona.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/db/builder/mona.py module
"""


import unittest
import os
import json
from tempfile import TemporaryDirectory
import contextlib
import io
import pathlib
import glob

from idpp.db.builder.mona import (
   mona_chunks_exist, chunk_mona_json, _parse_mona_entry, add_mona_chunks_to_idppdb
)
from idpp.db.util import create_db, IdPPdb
from idpp.test.__include import TEST_INCLUDE_DIR


# define file path to mock_MoNA-export-Experimental_Spectra.json
_MOCK_MONA_SPECTRA_JSON = os.path.join(TEST_INCLUDE_DIR, "mock_MoNA-export-Experimental_Spectra.json")


class Test_MockMonaSpectraExists(unittest.TestCase):
    """ make sure that the built-in mock_MoNA-export-Experimental_Spectra.json exists """

    def test_MMSE_file_exists(self):
        """ built-in mock_MoNA-export-Experimental_Spectra.json file exists """
        self.assertTrue(os.path.isfile(_MOCK_MONA_SPECTRA_JSON))


class TestMonaChunksExist(unittest.TestCase):
    """ tests for the mona_chunks_exist function """

    def test_MCE_no_chunk_dir(self):
        """ if chunk_dir doesn't exist, should return False """
        self.assertFalse(mona_chunks_exist("this directory does not exist"))

    def test_MCE_chunk_dir_exists_but_empty(self):
        """ if chunk_dir exists but is empty, should return False """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            self.assertFalse(mona_chunks_exist(tmp_dir))

    def test_MCE_chunk_dir_and_chunk_files_exist(self):
        """ if chunk_dir exists and has enough of the JSON files, should return True """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # add 101 JSON files
            for i in range(102):
                pathlib.Path.touch(os.path.join(tmp_dir, f"{i + 1:04d}.json"))
            self.assertTrue(mona_chunks_exist(tmp_dir))


class TestChunkMonaJson(unittest.TestCase):
    """ tests for the chunk_mona_json function """

    def test_CMJ_no_hmdb_metabolites_file(self):
        """ if the hmdb_metabolites file does not exist, should raise an error """
        with self.assertRaises(FileNotFoundError), contextlib.redirect_stdout(io.StringIO()) as _:
            chunk_mona_json("this file does not exist", "chunks")
    
    def test_CMJ_chunk_successfull(self):
        """ test that chunking works as expected """
        # temporarily redirect stdout to suppress the print messages
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                contextlib.redirect_stdout(io.StringIO()) as _:
            chunk_dir = os.path.join(tmp_dir, "chunks")
            os.mkdir(chunk_dir)
            chunk_mona_json(_MOCK_MONA_SPECTRA_JSON, chunk_dir)
            # this mock input should generate at least 4 chunks
            self.assertGreaterEqual(len(glob.glob("[0-9][0-9][0-9][0-9].json", 
                                                  root_dir=chunk_dir)), 2)
    
    def test_CMJ_chunk_dir_does_not_exist(self):
        """ if chunk_dir does not exist, it should be created then things work as normal """
        # temporarily redirect stdout to suppress the print messages
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                contextlib.redirect_stdout(io.StringIO()) as _:
            chunk_dir = os.path.join(tmp_dir, "chunks")
            #os.mkdir(chunk_dir)  <- let chunk_hmdb_xml create the directory
            chunk_mona_json(_MOCK_MONA_SPECTRA_JSON, chunk_dir)
            # this mock input should generate at least 4 chunks
            self.assertGreaterEqual(len(glob.glob("[0-9][0-9][0-9][0-9].json", 
                                                  root_dir=chunk_dir)), 2)


class Test_ParseMonaEntry(unittest.TestCase):
    """ tests for the _parse_mona_entry function """

    def test_PME_parsing_ok(self):
        """ validate _parse_mona_entry outputs using the mock MoNA export data """
        with open(_MOCK_MONA_SPECTRA_JSON, "r", encoding="utf8") as jf:
            entries = json.load(jf)
        valid_entries = 0
        for entry in entries:        
            if (parsed := _parse_mona_entry(entry)) is not None:
                # count up the number of valid parsed entries
                valid_entries += 1
                # make sure all of the required fields have been filled
                p_req, _ = parsed
                self.assertIsNotNone(p_req.names)
                self.assertGreater(len(p_req.names), 0)
                self.assertIsNotNone(p_req.mz)
                self.assertIsNotNone(p_req.adduct)
                self.assertIsNotNone(p_req.spectrum)
        # make sure enough valid entries were found out of this mock dataset
        # with 1000 entries total (grabbed from a chunk in the middle of the real file)
        # oh weird, it turns out all of them are parsable, lol ok that's the requirement then 
        self.assertGreaterEqual(valid_entries, 1000)


class TestAddMonaChunksToIdppdb(unittest.TestCase):
    """ tests for add_mona_chunks_to_idppdb """

    def test_AMCTI_no_chunk_dir(self):
        """ if the chunk directory with JSON files does not exist, should raise an error """
        with self.assertRaises(RuntimeError):
            with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                    contextlib.redirect_stdout(io.StringIO()) as _:
                dbf = os.path.join(tmp_dir, "idpp.db")
                create_db(dbf)
                db = IdPPdb(dbf)
                add_mona_chunks_to_idppdb(db, "chunk_dir_does_not_exist")

    def test_AMCTI_add_mock_data(self):
        """ add mock MoNA experimental MS/MS dataset to IdPPdb """
        # temporarily redirect stdout to suppress the print messages
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                contextlib.redirect_stdout(io.StringIO()) as _:
            # init the database
            dbf = os.path.join(tmp_dir, "idpp.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            # make some mock chunk data
            chunk_dir = os.path.join(tmp_dir, "chunks")
            os.mkdir(chunk_dir)
            chunk_mona_json(_MOCK_MONA_SPECTRA_JSON, chunk_dir)
            add_mona_chunks_to_idppdb(db, chunk_dir)
            # TODO: Check the counts of some of the tables?

    def test_AMCTI_add_mock_data_no_combine(self):
        """ add mock MoNA experimental MS/MS dataset to IdPPdb without combining spectra """
        # temporarily redirect stdout to suppress the print messages
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                contextlib.redirect_stdout(io.StringIO()) as _:
            # init the database
            dbf = os.path.join(tmp_dir, "idpp.db")
            create_db(dbf)
            db = IdPPdb(dbf, combine_ms2=False)
            # make some mock chunk data
            chunk_dir = os.path.join(tmp_dir, "chunks")
            os.mkdir(chunk_dir)
            chunk_mona_json(_MOCK_MONA_SPECTRA_JSON, chunk_dir)
            add_mona_chunks_to_idppdb(db, chunk_dir)
            # TODO: Check the counts of some of the tables?


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsMona = unittest.TestSuite()
AllTestsMona.addTests([
    _loader.loadTestsFromTestCase(Test_MockMonaSpectraExists),
    _loader.loadTestsFromTestCase(TestMonaChunksExist),
    _loader.loadTestsFromTestCase(TestChunkMonaJson),
    _loader.loadTestsFromTestCase(Test_ParseMonaEntry),
    _loader.loadTestsFromTestCase(TestAddMonaChunksToIdppdb),
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsMona)

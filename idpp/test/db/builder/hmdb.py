"""
    idpp/test/db/builder/hmdb.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/db/builder/hmdb.py module
"""


import os
import glob
import pathlib
import contextlib
import io
import unittest
from tempfile import TemporaryDirectory

from idpp.db.util import create_db, IdPPdb
from idpp.test.__include import TEST_INCLUDE_DIR
from idpp.db.builder.hmdb import (
    hmdb_chunks_exist, chunk_hmdb_xml, add_hmdb_chunks_to_idppdb
)


# define file path to mock_hmdb_metabolites.xml
_MOCK_HMDB_METABS_XML = os.path.join(TEST_INCLUDE_DIR, "mock_hmdb_metabolites.xml")


class Test_MockHmdbMetabolitesExists(unittest.TestCase):
    """ make sure that the built-in mock_hmdb_metabolites.xml exists """

    def test_MHME_file_exists(self):
        """ built-in mock_hmdb_metabolites.xml file exists """
        self.assertTrue(os.path.isfile(_MOCK_HMDB_METABS_XML))


class TestHmdbChunksExist(unittest.TestCase):
    """ tests for the hmdb_chunks_exist function """

    def test_HCE_no_chunk_dir(self):
        """ if chunk_dir doesn't exist, should return False """
        self.assertFalse(hmdb_chunks_exist("this directory does not exist"))

    def test_HCE_chunk_dir_exists_but_empty(self):
        """ if chunk_dir exists but is empty, should return False """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            self.assertFalse(hmdb_chunks_exist(tmp_dir))

    def test_HCE_chunk_dir_and_chunk_files_exist(self):
        """ if chunk_dir exists and has enough of the XML files, should return True """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # add 101 XML files
            for i in range(102):
                pathlib.Path.touch(os.path.join(tmp_dir, f"{i + 1:04d}.xml"))
            self.assertTrue(hmdb_chunks_exist(tmp_dir))


class TestChunkHmdbXml(unittest.TestCase):
    """ tests for the chunk_hmdb_xml function """

    def test_CHX_no_hmdb_metabolites_file(self):
        """ if the hmdb_metabolites file does not exist, should raise an error """
        with self.assertRaises(FileNotFoundError), contextlib.redirect_stdout(io.StringIO()) as _:
            chunk_hmdb_xml("this file does not exist", "chunks")
    
    def test_CHX_chunk_successfull(self):
        """ test that chunking works as expected """
        # temporarily redirect stdout to suppress the print messages
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                contextlib.redirect_stdout(io.StringIO()) as _:
            chunk_dir = os.path.join(tmp_dir, "chunks")
            os.mkdir(chunk_dir)
            chunk_hmdb_xml(_MOCK_HMDB_METABS_XML, chunk_dir)
            # this mock input should generate at least 4 chunks
            self.assertGreaterEqual(len(glob.glob("[0-9][0-9][0-9][0-9].xml", 
                                                  root_dir=chunk_dir)), 4)
    
    def test_CHX_chunk_dir_does_not_exist(self):
        """ if chunk_dir does not exist, it should be created then things work as normal """
        # temporarily redirect stdout to suppress the print messages
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                contextlib.redirect_stdout(io.StringIO()) as _:
            chunk_dir = os.path.join(tmp_dir, "chunks")
            #os.mkdir(chunk_dir)  <- let chunk_hmdb_xml create the directory
            chunk_hmdb_xml(_MOCK_HMDB_METABS_XML, chunk_dir)
            # this mock input should generate at least 4 chunks
            self.assertGreaterEqual(len(glob.glob("[0-9][0-9][0-9][0-9].xml", 
                                                  root_dir=chunk_dir)), 4)


class TestAddHmdbChunksToIdppdb(unittest.TestCase):
    """ tests for add_hmdb_chunks_to_idppdb """

    def test_AHCTI_no_chunk_dir(self):
        """ if the chunk directory does not exist, should raise an error """
        with self.assertRaises(RuntimeError):
            with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                    contextlib.redirect_stdout(io.StringIO()) as _:
                dbf = os.path.join(tmp_dir, "idpp.db")
                create_db(dbf)
                db = IdPPdb(dbf)
                add_hmdb_chunks_to_idppdb(db, "chunk dir does not exist")

    def test_AHCTI_add_mock_data(self):
        """ add mock HMDB dataset to IdPPdb """
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
            chunk_hmdb_xml(_MOCK_HMDB_METABS_XML, chunk_dir)
            add_hmdb_chunks_to_idppdb(db, chunk_dir)
            # TODO: Check the counts of some of the tables?


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsHmdb = unittest.TestSuite()
AllTestsHmdb.addTests([
    _loader.loadTestsFromTestCase(Test_MockHmdbMetabolitesExists),
    _loader.loadTestsFromTestCase(TestHmdbChunksExist),
    _loader.loadTestsFromTestCase(TestChunkHmdbXml),
    _loader.loadTestsFromTestCase(TestAddHmdbChunksToIdppdb)
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsHmdb)

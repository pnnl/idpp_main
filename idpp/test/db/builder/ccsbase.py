"""
    idpp/test/db/builder/ccsbase.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/db/builder/ccsbase.py module
"""


import unittest
import os
from tempfile import TemporaryDirectory
import contextlib
import io
import sqlite3

from idpp.test.__include import TEST_INCLUDE_DIR
from idpp.db.util import create_db, IdPPdb
from idpp.db.builder.ccsbase import (
    _ccsbase_iter, add_ccsbase_to_idppdb
)


# define file path to mock_C3S.db
_MOCK_C3SDB_FILE = os.path.join(TEST_INCLUDE_DIR, "mock_C3S.db")


class Test_MockC3sdbExists(unittest.TestCase):
    """ make sure that the built-in mock_C3S.db exists """

    def test_MCE_file_exists(self):
        """ built-in mock_C3S.db file exists """
        self.assertTrue(os.path.isfile(_MOCK_C3SDB_FILE))


class Test_CcsbaseIter(unittest.TestCase):
    """ tests for the _ccsbase_iter function """

    def test_CI_iter_mock_data(self):
        """ test iterating rows in mock data, should be no errors """
        con = sqlite3.connect(_MOCK_C3SDB_FILE)
        cur = con.cursor()
        for (g_id, name, adduct, z, mz, ccs, smi, 
             src_tag, ccs_type, ccs_method) in _ccsbase_iter(cur):
            pass
            # TODO: validate the yielded values?
        con.close()


# class Test_MakeSrc(unittest.TestCase):
#     """ tests for the _make_src function """
    

class TestAddCcsbaseToIdppdb(unittest.TestCase):
    """ tests for the add_ccsbase_to_idppdb function """

    def test_ACTI_c3sdb_file_not_found(self):
        """ should raise an error if the C3S.db is not found """
        with self.assertRaises(FileNotFoundError):
            with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                    contextlib.redirect_stdout(io.StringIO()) as _:
                dbf = os.path.join(tmp_dir, "idpp.db")
                create_db(dbf)
                db = IdPPdb(dbf)
                add_ccsbase_to_idppdb(db, "this file does not exist")
    
    def test_ACTI_add_mock_data(self):
        """ add mock CCSbase to IdPPdb """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                contextlib.redirect_stdout(io.StringIO()) as _:
            # init the database
            dbf = os.path.join(tmp_dir, "idpp.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            add_ccsbase_to_idppdb(db, _MOCK_C3SDB_FILE)
            # TODO: Check the counts of some of the tables?


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsCcsbase = unittest.TestSuite()
AllTestsCcsbase.addTests([
    _loader.loadTestsFromTestCase(Test_MockC3sdbExists),
    _loader.loadTestsFromTestCase(Test_CcsbaseIter),
    _loader.loadTestsFromTestCase(TestAddCcsbaseToIdppdb),
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsCcsbase)

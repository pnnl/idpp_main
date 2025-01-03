"""
    idpp/test/db/builder/nist20.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/db/builder/nist20.py module
"""


import os
import contextlib
import io
import unittest
from tempfile import TemporaryDirectory

from idpp.db.util import create_db, IdPPdb
from idpp.test.__include import TEST_INCLUDE_DIR
from idpp.db.builder.nist20 import (
    add_nist20_msms_to_idppdb
)


# define file path to mock_msms_2020.db
_MOCK_NIST20_DB = os.path.join(TEST_INCLUDE_DIR, "mock_msms_2020.db")


class Test_MockNist20DbExists(unittest.TestCase):
    """ make sure that the built-in mock_msms_2020.db exists """

    def test_MNDE_file_exists(self):
        """ built-in mock_msms_2020.db file exists """
        self.assertTrue(os.path.isfile(_MOCK_NIST20_DB))


class TestAddNist20MsmsToIdppdb(unittest.TestCase):
    """ tests for add_nist20_msms_to_idppdb """

    def test_ANMTI_no_db_file(self):
        """ if the database file does not exist, should raise an error """
        with self.assertRaises(FileNotFoundError):
            with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                    contextlib.redirect_stdout(io.StringIO()) as _:
                dbf = os.path.join(tmp_dir, "idpp.db")
                create_db(dbf)
                db = IdPPdb(dbf)
                add_nist20_msms_to_idppdb(db, "the database does not exist")

    def test_ANMTI_add_mock_data(self):
        """ add mock MoNA experimental MS/MS dataset to IdPPdb """
        # temporarily redirect stdout to suppress the print messages
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                contextlib.redirect_stdout(io.StringIO()) as _:
            # init the database
            dbf = os.path.join(tmp_dir, "idpp.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            add_nist20_msms_to_idppdb(db, _MOCK_NIST20_DB)
            # TODO: Check the counts of some of the tables?

    def test_ANMTI_add_mock_data_no_combine(self):
        """ add mock MoNA experimental MS/MS dataset to IdPPdb without combining spectra """
        # temporarily redirect stdout to suppress the print messages
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                contextlib.redirect_stdout(io.StringIO()) as _:
            # init the database
            dbf = os.path.join(tmp_dir, "idpp.db")
            create_db(dbf)
            db = IdPPdb(dbf, combine_ms2=False)
            add_nist20_msms_to_idppdb(db, _MOCK_NIST20_DB)
            # TODO: Check the counts of some of the tables?


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsNist20 = unittest.TestSuite()
AllTestsNist20.addTests([
    _loader.loadTestsFromTestCase(Test_MockNist20DbExists),
    _loader.loadTestsFromTestCase(TestAddNist20MsmsToIdppdb)
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsNist20)

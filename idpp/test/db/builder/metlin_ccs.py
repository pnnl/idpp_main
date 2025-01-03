"""
    idpp/test/db/builder/ccs_compendium.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/db/builder/ccs_compendium.py module
"""


import unittest
import os
from tempfile import TemporaryDirectory
import contextlib
import io

from idpp.test.__include import TEST_INCLUDE_DIR
from idpp.db.util import create_db, IdPPdb
from idpp.db.builder.metlin_ccs import (
    _iter_metlin_ccs,
    add_metlin_ccs_to_idppdb
)


# define file path to mock_?
_MOCK_METLIN_FILE = os.path.join(TEST_INCLUDE_DIR, "mock_METLIN-CCS-03-15-2024.xlsx")


class Test_MockMetlinFileExists(unittest.TestCase):
    """ make sure that the built-in mock_METLIN-CCS-03-15-2024.xlsx exists """

    def test_MMFE_file_exists(self):
        """ built-in mock_METLIN-CCS-03-15-2024.xlsx file exists """
        self.assertTrue(os.path.isfile(_MOCK_METLIN_FILE))


class Test_IterMetlinCcs(unittest.TestCase):
    """ tests for the _iter_metlin_ccs function """

    def test_IMC_iter_mock_data(self):
        """ test iterating rows in mock data, should be no errors """
        for (name, 
             adduct, 
             formula, 
             metlin_id, 
             mz, 
             z,
             ccs) in _iter_metlin_ccs(_MOCK_METLIN_FILE, (0.2692, 121.54)):
            pass
            # TODO: validate the yielded values?


class TestAddMetlinCcsToIdppdb(unittest.TestCase):
    """ tests for the add_metlin_ccs_to_idppdb function """

    def test_AMCTI_mock_file_not_found(self):
        """ should raise an error if the mock data file is not found """
        with self.assertRaises(FileNotFoundError):
            with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                    contextlib.redirect_stdout(io.StringIO()) as _:
                dbf = os.path.join(tmp_dir, "idpp.db")
                create_db(dbf)
                db = IdPPdb(dbf)
                add_metlin_ccs_to_idppdb(db, "this file does not exist")
    
    def test_AMCTI_add_mock_data(self):
        """ add mock METLIN-CCS dataset to IdPPdb """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                contextlib.redirect_stdout(io.StringIO()) as _:
            # init the database
            dbf = os.path.join(tmp_dir, "idpp.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            add_metlin_ccs_to_idppdb(db, _MOCK_METLIN_FILE)
            # TODO: Check the counts of some of the tables?


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsMetlinCcs = unittest.TestSuite()
AllTestsMetlinCcs.addTests([
    _loader.loadTestsFromTestCase(Test_MockMetlinFileExists),
    _loader.loadTestsFromTestCase(Test_IterMetlinCcs),
    _loader.loadTestsFromTestCase(TestAddMetlinCcsToIdppdb),
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsMetlinCcs)

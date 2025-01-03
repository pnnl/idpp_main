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
from idpp.db.builder.ccs_compendium import (
    _compendium_iter, add_ccs_compendium_to_idppdb
)


# define file path to mock_UnifiedCCSCompendium
_MOCK_COMPENDIUM_FILE = os.path.join(TEST_INCLUDE_DIR, "mock_UnifiedCCSCompendium.csv")


class Test_MockCcsCompendiumExists(unittest.TestCase):
    """ make sure that the built-in mock_UnifiedCCSCompendium exists """

    def test_MCCE_file_exists(self):
        """ built-in mock_UnifiedCCSCompendium file exists """
        self.assertTrue(os.path.isfile(_MOCK_COMPENDIUM_FILE))


class Test_CompendiumIter(unittest.TestCase):
    """ tests for the _compendium_iter function """

    def test_iter_mock_data(self):
        """ test iterating rows in mock data, should be no errors """
        for (name, form, inchi, inchikey, 
                cls_info, mz, adduct, z, ccs) in _compendium_iter(_MOCK_COMPENDIUM_FILE):
            pass
            # TODO: validate the yielded values?



class TestAddCcsCompendiumToIdppdb(unittest.TestCase):
    """ tests for the add_ccs_compendium_to_idppdb function """

    def test_ACCTI_rtdata_file_not_found(self):
        """ should raise an error if the raw_data directory is not found """
        with self.assertRaises(FileNotFoundError):
            with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                    contextlib.redirect_stdout(io.StringIO()) as _:
                dbf = os.path.join(tmp_dir, "idpp.db")
                create_db(dbf)
                db = IdPPdb(dbf)
                add_ccs_compendium_to_idppdb(db, "this file does not exist")
    
    def test_ACCTI_add_mock_data(self):
        """ add mock CCS compendium to IdPPdb """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir, \
                contextlib.redirect_stdout(io.StringIO()) as _:
            # init the database
            dbf = os.path.join(tmp_dir, "idpp.db")
            create_db(dbf)
            db = IdPPdb(dbf)
            add_ccs_compendium_to_idppdb(db, _MOCK_COMPENDIUM_FILE)
            # TODO: Check the counts of some of the tables?


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsCcsCompendium = unittest.TestSuite()
AllTestsCcsCompendium.addTests([
    _loader.loadTestsFromTestCase(Test_MockCcsCompendiumExists),
    _loader.loadTestsFromTestCase(Test_CompendiumIter),
    _loader.loadTestsFromTestCase(TestAddCcsCompendiumToIdppdb),
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsCcsCompendium)

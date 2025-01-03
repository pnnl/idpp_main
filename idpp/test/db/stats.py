"""
    idpp/test/db/stats.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/db/stats.py module
"""


import unittest

from idpp.db.stats import (
    compound_property_coverage, 
    property_source_distributions, 
    _main
)


class TestCompoundPropertyCoverage(unittest.TestCase):
    """ tests for the adduct_property_coverage function """

    def test_NO_TESTS_IMPLEMENTED_YET(self):
        """ placeholder, remove this function and implement tests """
        assert False, "not implemented yet"


class TestPropertySourceDistributions(unittest.TestCase):
    """ tests for the property_source_distributions function """

    def test_NO_TESTS_IMPLEMENTED_YET(self):
        """ placeholder, remove this function and implement tests """
        assert False, "not implemented yet"


class Test_Main(unittest.TestCase):
    """ tests for the _main function """
    
    def test_NO_TESTS_IMPLEMENTED_YET(self):
        """ placeholder, remove this function and implement tests """
        # TODO: mock CLI args that would normally get processed
        #       in _main when idpp/db/stats.py module is invoked directly
        assert False, "not implemented yet"


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsStats = unittest.TestSuite()
AllTestsStats.addTests([
    _loader.loadTestsFromTestCase(TestCompoundPropertyCoverage),
    _loader.loadTestsFromTestCase(TestPropertySourceDistributions),
    _loader.loadTestsFromTestCase(Test_Main)
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsStats)

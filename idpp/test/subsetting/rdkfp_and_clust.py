"""
    idpp/test/subsetting/rdkfp_and_clust.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/subsetting/rdkfp_and_clust.py module
"""


import unittest

from idpp.subsetting import (
    rdkfp_and_clust
)


class TestRdkfpAndClust(unittest.TestCase):
    """ tests for the rdkfp_and_clust function """

    def test_NO_TESTS_IMPLEMENTED_YET(self):
        """ placeholder, remove this function and implement tests """
        assert False, "not implemented yet"


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsRdkfpAndClust = unittest.TestSuite()
AllTestsRdkfpAndClust.addTests([
    _loader.loadTestsFromTestCase(TestRdkfpAndClust)
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsRdkfpAndClust)


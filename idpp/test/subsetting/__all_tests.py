"""
    idpp/test/subsetting/__all_tests.py

    Dylan Ross (dylan.ross@pnnl.gov)

    special module for grouping all TestSuites from submodules/subpackages
"""


import unittest

from idpp.test.subsetting.rdkfp_and_clust import AllTestsRdkfpAndClust


# collect tests from this subpackage
AllTestsSubsetting = unittest.TestSuite()
AllTestsSubsetting.addTests([
    AllTestsRdkfpAndClust
])

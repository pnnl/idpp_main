"""
    idpp/test/probability/__all_tests.py

    Dylan Ross (dylan.ross@pnnl.gov)

    special module for grouping all TestSuites from submodules/subpackages
"""


import unittest

from idpp.test.probability.trees import AllTestsTrees
from idpp.test.probability.analysis import AllTestsAnalysis


# collect tests from this subpackage
AllTestsProbability = unittest.TestSuite()
AllTestsProbability.addTests([
    AllTestsTrees,
    AllTestsAnalysis
])


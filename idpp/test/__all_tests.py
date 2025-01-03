"""
    idpp/test/__all_tests.py

    Dylan Ross (dylan.ross@pnnl.gov)

    special module for grouping all TestSuites from submodules/subpackages
"""


import unittest

from idpp.test.db.__all_tests import AllTestsDb
from idpp.test.ions import AllTestsIons
from idpp.test.msms.__all_tests import AllTestsMsms
from idpp.test.probability.__all_tests import AllTestsProbability
from idpp.test.subsetting.__all_tests import AllTestsSubsetting


# collect tests
AllTests = unittest.TestSuite()
AllTests.addTests([
    AllTestsDb,
    AllTestsIons,
    AllTestsMsms,
    AllTestsProbability,
    AllTestsSubsetting
])

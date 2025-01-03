"""
    idpp/test/db/__all_tests.py

    Dylan Ross (dylan.ross@pnnl.gov)

    special module for grouping all TestSuites from submodules/subpackages
"""


import unittest

from idpp.test.db.builder.__all_tests import AllTestsBuilder
from idpp.test.db.stats import AllTestsStats
from idpp.test.db.util import AllTestsUtil


# collect tests from this subpackage
AllTestsDb = unittest.TestSuite()
AllTestsDb.addTests([
    AllTestsBuilder,
    AllTestsStats,
    AllTestsUtil
])


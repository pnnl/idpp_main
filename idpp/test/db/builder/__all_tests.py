"""
    idpp/test/db/builder/__all_tests.py

    Dylan Ross (dylan.ross@pnnl.gov)

    special module for grouping all TestSuites from submodules/subpackages
"""


import unittest

from idpp.test.db.builder._util import AllTests_Util
from idpp.test.db.builder.ccs_compendium import AllTestsCcsCompendium
from idpp.test.db.builder.ccsbase import AllTestsCcsbase
from idpp.test.db.builder.hmdb import AllTestsHmdb
from idpp.test.db.builder.mona import AllTestsMona
from idpp.test.db.builder.nist20 import AllTestsNist20
from idpp.test.db.builder.report import AllTestsReport
from idpp.test.db.builder.metlin_ccs import AllTestsMetlinCcs


# collect tests from this subpackage
AllTestsBuilder = unittest.TestSuite()
AllTestsBuilder.addTests([
    AllTests_Util,
    AllTestsCcsCompendium,
    AllTestsCcsbase,
    AllTestsHmdb,
    AllTestsMona,
    AllTestsNist20,
    AllTestsReport,
    AllTestsMetlinCcs
])


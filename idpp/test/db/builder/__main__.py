"""
    idpp/test/db/builder/__main__.py

    Dylan Ross (dylan.ross@pnnl.gov)

    runs unit tests for builder subpackage
"""


import unittest

from idpp.test.db.builder.__all_tests import AllTestsBuilder


# run all defined TestCases for this subpackage
runner = unittest.TextTestRunner(verbosity=2)
result = runner.run(AllTestsBuilder)

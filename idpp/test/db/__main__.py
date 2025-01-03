"""
    idpp/test/db/__main__.py

    Dylan Ross (dylan.ross@pnnl.gov)

    runs unit tests for db subpackage
"""


import unittest

from idpp.test.db.__all_tests import AllTestsDb


# run all defined TestCases for this subpackage
runner = unittest.TextTestRunner(verbosity=2)
result = runner.run(AllTestsDb)

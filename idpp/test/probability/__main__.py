"""
    idpp/test/probability/__main__.py

    Dylan Ross (dylan.ross@pnnl.gov)

    runs unit tests for probability subpackage
"""


import unittest

from idpp.test.probability.__all_tests import AllTestsProbability


# run all defined TestCases for this subpackage
runner = unittest.TextTestRunner(verbosity=2)
result = runner.run(AllTestsProbability)

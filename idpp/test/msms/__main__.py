"""
    idpp/test/msms/__main__.py

    Dylan Ross (dylan.ross@pnnl.gov)

    runs unit tests for msms subpackage
"""


import unittest

from idpp.test.msms.__all_tests import AllTestsMsms


# run all defined TestCases for this subpackage
runner = unittest.TextTestRunner(verbosity=2)
result = runner.run(AllTestsMsms)

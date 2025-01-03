"""
    idpp/test/subsetting/__main__.py

    Dylan Ross (dylan.ross@pnnl.gov)

    runs unit tests for subsetting subpackage
"""


import unittest

from idpp.test.subsetting import AllTestsSubsetting


# run all defined TestCases for this subpackage
runner = unittest.TextTestRunner(verbosity=2)
result = runner.run(AllTestsSubsetting)

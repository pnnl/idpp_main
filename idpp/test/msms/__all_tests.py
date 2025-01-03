"""
    idpp/test/msms/__all_tests.py

    Dylan Ross (dylan.ross@pnnl.gov)

    special module for grouping all TestSuites from submodules/subpackages
"""


import unittest

from idpp.test.msms.spectra import AllTestsSpectra


# collect tests from this subpackage
AllTestsMsms = unittest.TestSuite()
AllTestsMsms.addTests([
    AllTestsSpectra
])




"""
    idpp/test/probability/analysis.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/probability/analysis.py module
"""


import unittest


class TestPlaceholder(unittest.TestCase):
    """ tests for the ? """

    def test_placeholder(self):
        assert False, "none implemented yet"


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsAnalysis = unittest.TestSuite()
AllTestsAnalysis.addTests([
    _loader.loadTestsFromTestCase(TestPlaceholder)
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsAnalysis)

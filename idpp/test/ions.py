"""
    idpp/test/ions.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/ions.py module
"""


import unittest

from idpp.ions import (
    _get_replacements, _protonate, _ion_adduct, ionize_smi
)


class Test_GetReplacements(unittest.TestCase):
    """ tests for the _get_replacements function """

    def test_NO_TESTS_IMPLEMENTED_YET(self):
        """ placeholder, remove this function and implement tests """
        assert False, "not implemented yet"


class Test_Protonate(unittest.TestCase):
    """ tests for the _protonate function """

    def test_NO_TESTS_IMPLEMENTED_YET(self):
        """ placeholder, remove this function and implement tests """
        assert False, "not implemented yet"


class Test_IonAdduct(unittest.TestCase):
    """ tests for the _ion_adduct function """

    def test_NO_TESTS_IMPLEMENTED_YET(self):
        """ placeholder, remove this function and implement tests """
        assert False, "not implemented yet"


class TestIonizeSmi(unittest.TestCase):
    """ tests for the ionize_smi function """

    def test_NO_TESTS_IMPLEMENTED_YET(self):
        """ placeholder, remove this function and implement tests """
        assert False, "not implemented yet"


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsIons = unittest.TestSuite()
AllTestsIons.addTests([
    _loader.loadTestsFromTestCase(Test_GetReplacements),
    _loader.loadTestsFromTestCase(Test_Protonate),
    _loader.loadTestsFromTestCase(Test_IonAdduct),
    _loader.loadTestsFromTestCase(TestIonizeSmi)
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsIons)


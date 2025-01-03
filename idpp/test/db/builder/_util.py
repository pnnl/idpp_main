"""
    idpp/test/db/builder/_util.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/db/builder/_util.py module
"""


import unittest

import numpy as np

from idpp.db.builder._util import (
    parse_ce, str_to_ms2
)


class TestParseCe(unittest.TestCase):
    """ tests for the parse_ce function """

    def test_PS_ce_str_is_none(self):
        """ when None is passed in as ce_str, should return None """
        self.assertIsNone(parse_ce(None))

    def test_PC_unparsable(self):
        """ test some cases that should not be parsable and should return None """
        for ce_str in [
            "not parsable", 
            "normalized 34V", 
            "50%", 
            "50 %"
            "",
        ]:
            self.assertIsNone(parse_ce(ce_str), f"{ce_str} did not return None")

    def test_PC_parsable(self):
        """ test some cases that should be parsable and should return correct value """
        for pre in ["", "CE", "CE-", "CE ", "NCE", "NCE-", "NCE "]:
            for v, expected in [("1", 1),
                      ("1.", 1),
                      ("1.0", 1),
                      ("1.1", 1),
                      ("69", 69),
                      ("69.", 69),
                      ("69.0", 69),
                      ("69.1", 69),
                      ("420", 420),
                      ("420.", 420),
                      ("420.0", 420),
                      ("420.1", 420)]:
                for post in ["", "V", "eV", " V", " eV", "NCE", " NCE"]:
                    ce_str = pre + v + post
                    self.assertEqual(parse_ce(ce_str), expected,
                                     f"{ce_str} did not return expected value (69)")


class TestStrToMS2(unittest.TestCase):
    """ tests for the str_to_ms2 function """

    def test_STM_bad_str_formats(self):
        """ bad spectrum string formats should raise ValueErrors """
        bads = [
            "", "bad spectrum string",
            "1", "1:", "1.:", 
            "1.:1 1", "1.:1 1:", "1.:1 1.:"
        ]
        for bad in bads:
            with self.assertRaises(ValueError,
                                   msg=f"spectrum string '{bad}' should have cause ValueError"):
                mzs, iis = str_to_ms2(bad)
    
    def _arrays_length_and_content_match(self, a, a_exp, a_label):
        """
        helper method that compares the lengths and contents of a and a_exp
        """
        self.assertEqual(len(a), len(a_exp),
                         msg=f"incorrect length of {a_label} array")
        for v, v_exp in zip(a, a_exp):
            self.assertEqual(v, v_exp,
                             msg=f"element {v} of {a_label} array does not match expected value ({v_exp})")

    def test_STM_expected_outputs(self):
        """ good spectrum strings should produce the expected output arrays """
        expected = [
            # (spectrum string, output arrays)
            ("1:1", (np.array([1.]), np.array([1.]))),
            ("1.:1", (np.array([1.]), np.array([1.]))),
            ("1.:1.", (np.array([1.]), np.array([1.]))),
            ("1e1:1e1", (np.array([10.]), np.array([10.]))),
            ("1.e1:1.e1", (np.array([10.]), np.array([10.]))),
            ("1.:1.0", (np.array([1.]), np.array([1.]))),
            ("1.:1 2.:2", (np.array([1., 2.]), np.array([1., 2.]))),
            ("1.:1   2.:2", (np.array([1., 2.]), np.array([1., 2.]))),
            ("1.1:1 2.2:2", (np.array([1.1, 2.2]), np.array([1., 2.]))),
            ("1.11111:1   2.22222:2", (np.array([1.11111, 2.22222]), np.array([1., 2.]))),
        ]
        for s, (exp_mzs, exp_iis) in expected:
            mzs, iis = str_to_ms2(s)
            self._arrays_length_and_content_match(mzs, exp_mzs, "mzs")
            self._arrays_length_and_content_match(iis, exp_iis, "iis")


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTests_Util = unittest.TestSuite()
AllTests_Util.addTests([
    _loader.loadTestsFromTestCase(TestParseCe),
    _loader.loadTestsFromTestCase(TestStrToMS2)
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTests_Util)

"""
    idpp/test/msms/spectra.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/msms/spectra.py module
"""


import unittest
import cProfile
import pstats

import numpy as np

from idpp.msms.spectra import (
    spec_norm, spec_entropy, spec_combine, spec_entropy_similarity
)


class TestSpecNorm(unittest.TestCase):
    """ tests for the spec_norm function """

    def test_NO_TESTS_IMPLEMENTED_YET(self):
        """ placeholder, remove this function and implement tests """
        assert False, "not implemented yet"


class TestSpecEntropy(unittest.TestCase):
    """ tests for the spec_entropy function """

    def test_NO_TESTS_IMPLEMENTED_YET(self):
        """ placeholder, remove this function and implement tests """
        assert False, "not implemented yet"


def _old_spec_combine(spectra, weights, comb_mztol=0.05):
    """ 
    old implementation of spec_combine
    """
    comb_mz, comb_i = [], []
    # iterate through spectra in set
    # combine m/zs that are within tolerance of one another
    # otherwise just add the new points
    for spectrum, weight in zip(spectra, weights):
        for mz, i in zip(*spectrum):
            iw = i * weight
            added = False
            for j in range(len(comb_mz)):
                if abs(mz - comb_mz[j]) <= comb_mztol:
                    comb_i[j] += iw
                    added = True
            if not added:
                comb_mz.append(mz)
                comb_i.append(iw)
    # sort by m/z
    idx = np.argsort(comb_mz)
    # return spectrum sorted by m/z
    # do not normalize
    return np.array(comb_mz)[idx], np.array(comb_i)[idx]


class TestSpecCombine(unittest.TestCase):
    """ tests for the spec_combine function """

    def _gen_spectrum(self):
        """ helper method for generating a mock MS/MS spectrum """
        # it can have anywhere from 20 to 100 peaks
        n_peaks = np.random.randint(20, 101)
        mz = np.random.choice(np.arange(50, 1000, 0.05), n_peaks)
        i = 100. * np.random.random(n_peaks)
        idx = np.argsort(mz)
        return mz[idx], i[idx]

    def _array_elements_equal(self, a, a_exp, test_label):
        """ make sure the elements of array a match array a_exp """
        for i, (v, v_exp) in enumerate(zip(a, a_exp)):
            self.assertEqual(v, v_exp,
                             msg=f"test: {test_label}: element {i} of array ({v}) does not match expected value ({v_exp})")

    def test_SC_equal_weights(self):
        """ combine spectra with equal weights """
        expected = [
            # (
            #     "label"
            #     [
            #         input_spectrum, 
            #         ...
            #     ], 
            #     [weights, ...], 
            #     combined_spectrum
            # ),
            (
                "two spectra with 3 peaks, 2 are overlapping",
                [
                    np.array([[123.456, 234.567, 345.678], [100., 200., 300.]]), 
                    np.array([[234.567, 345.678, 456.7890], [300., 200., 100.]])
                ], 
                [1., 1.], 
                np.array([[123.456, 234.567, 345.678, 456.789], [100, 500, 500, 100]])
            ),
            (
                "one empty spectrum, the other with 3 peaks",
                [
                    np.array([[], []]), 
                    np.array([[234.567, 345.678, 456.7890], [300., 200., 100.]])
                ], 
                [1., 1.], 
                np.array([[234.567, 345.678, 456.7890], [300., 200., 100.]])
            ),
        ]
        for test_label, input_spectra, weights, output_spectrum in expected:
            combined = spec_combine(input_spectra, weights)
            self._array_elements_equal(combined[0], output_spectrum[0], test_label)
            self._array_elements_equal(combined[1], output_spectrum[1], test_label)

    def test_SC_unequal_weights(self):
        """ combine spectra with unequal weights """
        expected = [
            # (
            #     "label"
            #     [
            #         input_spectrum, 
            #         ...
            #     ], 
            #     [weights, ...], 
            #     combined_spectrum
            # ),
            (
                "two spectra with 3 peaks, 2 are overlapping",
                [
                    np.array([[123.456, 234.567, 345.678], [100., 200., 300.]]), 
                    np.array([[234.567, 345.678, 456.7890], [300., 200., 100.]])
                ], 
                [1., 3.], 
                np.array([[123.456, 234.567, 345.678, 456.789], [100., 1100., 900., 300.]])
            ),
            (
                "one empty spectrum, the other with 3 peaks",
                [
                    np.array([[], []]), 
                    np.array([[234.567, 345.678, 456.7890], [300., 200., 100.]])
                ], 
                [1., 3.], 
                np.array([[234.567, 345.678, 456.7890], [900., 600., 300.]])
            ),
            (
                "two spectra with 3 peaks, 2 are overlapping",
                [
                    np.array([[234.567, 345.678, 456.7890], [300., 200., 100.]]),
                    np.array([[123.456, 234.567, 345.678], [100., 200., 300.]]) 
                ], 
                [3., 1.], 
                np.array([[123.456, 234.567, 345.678, 456.789], [100., 1100., 900., 300.]])
            ),
        ]
        for test_label, input_spectra, weights, output_spectrum in expected:
            combined = spec_combine(input_spectra, weights)
            self._array_elements_equal(combined[0], output_spectrum[0], test_label)
            self._array_elements_equal(combined[1], output_spectrum[1], test_label)

    def _test_SC_benchmark(self):
        """ benchmark performance of new and old implementations of spec_combine """
        print()
        np.random.seed(69420)
        n_spectra = 1000
        spectra = [self._gen_spectrum() for _ in range(n_spectra)]
        block_size = 100
        print("-" * 40)
        print("NEW")
        with cProfile.Profile() as pr:
            for block_idx in range(n_spectra // block_size):
                combined = spectra[block_idx * block_size]
                for other_spectrum in spectra[block_size * block_idx + 1:block_size * (block_idx + 1)]:
                    combined = spec_combine([combined, other_spectrum], [3., 3.])
            # report the stats
            stats = pstats.Stats(pr)
            stats.sort_stats(pstats.SortKey.CUMULATIVE)
            stats.print_stats()
        print("-" * 40)
        print("OLD")
        with cProfile.Profile() as pr:
            for block_idx in range(n_spectra // block_size):
                combined = spectra[block_idx * block_size]
                for other_spectrum in spectra[block_size * block_idx + 1:block_size * (block_idx + 1)]:
                    combined = _old_spec_combine([combined, other_spectrum], [3., 3.])
            # report the stats
            stats = pstats.Stats(pr)
            stats.sort_stats(pstats.SortKey.CUMULATIVE)
            stats.print_stats()


class TestSpecEntropySimilarity(unittest.TestCase):
    """ tests for the spec_entropy_similarity function """

    def test_NO_TESTS_IMPLEMENTED_YET(self):
        """ placeholder, remove this function and implement tests """
        assert False, "not implemented yet"


# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsSpectra = unittest.TestSuite()
AllTestsSpectra.addTests([
    _loader.loadTestsFromTestCase(TestSpecNorm),
    _loader.loadTestsFromTestCase(TestSpecEntropy),
    _loader.loadTestsFromTestCase(TestSpecCombine),
    _loader.loadTestsFromTestCase(TestSpecEntropySimilarity),
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsSpectra)


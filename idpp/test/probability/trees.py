"""
    idpp/test/probability/trees.py

    Dylan Ross (dylan.ross@pnnl.gov)

    tests for idpp/probability/trees.py module
"""


import unittest
from tempfile import TemporaryDirectory
import os

import numpy as np
from numpy import typing as npt

from idpp.probability.trees import (
    MzTree,
    CcsTree, 
    RtTree,
    load_tree
)


def _random_mzs(n: int) -> npt.NDArray[np.float64] :
    """ create a random array of m/z values of length n """
    return np.random.uniform(50., 1500., n)


def _random_ccss(n: int) -> npt.NDArray[np.float64] :
    """ create a random array of  CCSs of length n """
    return np.random.uniform(100., 500., n)


def _random_rts(n: int) -> npt.NDArray[np.float64] :
    """ create a random array of  RTs of length n """
    return np.random.uniform(0., 30., n)


class TestMzTree(unittest.TestCase):
    """ tests for the MzTree class """

    def test_MT_setup_and_query_noerr(self):
        """ set up MzTree with different sizes, query them, there should be no errors """
        for n in [1024, 2048, 4096, 8192, 16384]:
            mzt = MzTree(_random_mzs(n))
            _ = mzt.query_radius(10.)

    def test_MT_expected_query_radius_results(self):
        """ check that queries with different tolerances return the expected results """
        # set up some m/zs
        mzs = np.array([400., 400.001, 400.005, 400.01, 400.05, 400.1])
        mzt = MzTree(mzs, leaf_size=1)
        expected = [
            [
                [0], 
                [1], 
                [2], 
                [3], 
                [4], 
                [5]
            ],
            [
                [0, 1], 
                [0, 1], 
                [2], 
                [3], 
                [4], 
                [5]
            ],
            [
                [0, 1], 
                [0, 1, 2], 
                [1, 2], 
                [3], 
                [4], 
                [5]
            ],
            [
                [0, 1, 2, 3], 
                [0, 1, 2, 3], 
                [0, 1, 2, 3], 
                [0, 1, 2, 3], 
                [4], 
                [5]
            ],
        ]
        for exp, ppm in zip(expected, [1, 5, 10, 50]):
            res = mzt.query_radius(ppm)
            for i in range(6):
                self.assertListEqual(exp[i], res[i].tolist(),
                                     f"target index={i} does not match expected")
                
    def test_MT_save(self):
        """ test saving MzTree instance to file """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # init tree
            mzs = _random_mzs(1000)
            mzt = MzTree(mzs, leaf_size=1)
            # save the tree to file, should be no errors
            mzt.save(tmp_dir, 420)
            # make sure the file exists
            expected = os.path.join(tmp_dir, "MzTree_dsid=420.pkl")
            self.assertTrue(os.path.isfile(expected), 
                            f"expected saved tree file {expected} not found") 


class TestCcsTree(unittest.TestCase):
    """ tests for the CcsTree class """

    def test_CT_setup_and_query_noerr(self):
        """ set up CcsTree with different sizes, query them, there should be no errors """
        for n in [1024, 2048, 4096, 8192, 16384]:
            ccst = CcsTree(_random_ccss(n))
            _ = ccst.query_radius(10.)

    def test_CT_expected_query_radius_results(self):
        """ check that queries with different tolerances return the expected results """
        # set up some ccss
        ccss = np.array([200., 202., None, 204., 206., None])
        ccst = CcsTree(ccss, leaf_size=1)
        expected = [
            [
                [0, 1], 
                [0, 1, 3], 
                [], 
                [1, 3, 4], 
                [3, 4], 
                []
            ],
            [
                [0, 1, 3], 
                [0, 1, 3, 4], 
                [], 
                [0, 1, 3, 4], 
                [1, 3, 4], 
                []
            ],
            [
                [0, 1, 3, 4], 
                [0, 1, 3, 4], 
                [], 
                [0, 1, 3, 4], 
                [0, 1, 3, 4], 
                []
            ],
        ]
        for exp, percent in zip(expected, [1, 2, 3]):
            res = ccst.query_radius(percent)
            for i in range(6):
                self.assertListEqual(exp[i], res[i].tolist(),
                                     f"target index={i} does not match expected")

    def test_CT_save(self):
        """ test saving CcsTree instance to file """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # init tree
            ccss = _random_ccss(1000)
            ccst = CcsTree(ccss, leaf_size=1)
            # save the tree to file, should be no errors
            ccst.save(tmp_dir, 420)
            # make sure the file exists
            expected = os.path.join(tmp_dir, "CcsTree_dsid=420.pkl")
            self.assertTrue(os.path.isfile(expected), 
                            f"expected saved tree file {expected} not found")    


class TestRtTree(unittest.TestCase):
    """ tests for the RtTree class """

    def test_RT_setup_and_query_noerr(self):
        """ set up RtTree with different sizes, query them, there should be no errors """
        for n in [1024, 2048, 4096, 8192, 16384]:
            rtt = RtTree(_random_rts(n))
            _ = rtt.query_radius(10.)

    def test_RT_expected_query_radius_results(self):
        """ check that queries with different tolerances return the expected results """
        # set up some RTs
        rts = np.array([10., 10.1, 10.5, 11., 15., 20.])
        rtt = RtTree(rts, leaf_size=1)
        expected = [
            [
                [0], 
                [1], 
                [2], 
                [3], 
                [4], 
                [5]
            ],
            [
                [0, 1], 
                [0, 1], 
                [2], 
                [3], 
                [4], 
                [5]
            ],
            [
                [0, 1, 2], 
                [0, 1, 2], 
                [0, 1, 2, 3], 
                [2, 3], 
                [4], 
                [5]
            ],
            [
                [0, 1, 2, 3], 
                [0, 1, 2, 3], 
                [0, 1, 2, 3], 
                [0, 1, 2, 3], 
                [4], 
                [5]
            ],
        ]
        for exp, tol in zip(expected, [0.05, 0.1, 0.5, 1.]):
            res = rtt.query_radius(tol)
            for i in range(6):
                self.assertListEqual(exp[i], res[i].tolist(),
                                     f"target index={i} does not match expected")
                
    def test_RT_save(self):
        """ test saving RtTree instance to file """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # init tree
            rts = _random_rts(1000)
            rtt = RtTree(rts, leaf_size=1)
            # save the tree to file, should be no errors
            rtt.save(tmp_dir, 420)
            # make sure the file exists
            expected = os.path.join(tmp_dir, "RtTree_dsid=420.pkl")
            self.assertTrue(os.path.isfile(expected), 
                            f"expected saved tree file {expected} not found")


class TestMs2Tree(unittest.TestCase):
    """ tests for the Ms2Tree class """

    def test_M2T_placeholder(self):
        """ placeholder, remove this function and implement tests """
        assert False, "not implemented yet"


class TestLoadTree(unittest.TestCase):
    """ tests for the load_tree function """

    def test_LT_load_the_trees(self):
        """ create some trees with random data and make sure they can be reloaded successfully """
        with TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
            # init trees
            mzs = _random_mzs(1000)
            mzt = MzTree(mzs, leaf_size=1)
            ccss = _random_ccss(1000)
            ccst = CcsTree(ccss, leaf_size=1)
            rts = _random_rts(1000)
            rtt = RtTree(rts, leaf_size=1)
            # save the trees to file
            mzt.save(tmp_dir, 420)
            ccst.save(tmp_dir, 420)
            rtt.save(tmp_dir, 420)
            # load the trees, should be no errors
            mzt_2 = load_tree(os.path.join(tmp_dir, "MzTree_dsid=420.pkl"))
            rtt_2 = load_tree(os.path.join(tmp_dir, "RtTree_dsid=420.pkl"))
            ccst_2 = load_tree(os.path.join(tmp_dir, "CcsTree_dsid=420.pkl"))
    

# group all of the tests from this module into a TestSuite
_loader = unittest.TestLoader()
AllTestsTrees = unittest.TestSuite()
AllTestsTrees.addTests([
    _loader.loadTestsFromTestCase(TestMzTree),
    _loader.loadTestsFromTestCase(TestCcsTree),
    _loader.loadTestsFromTestCase(TestRtTree),
    _loader.loadTestsFromTestCase(TestMs2Tree),
    _loader.loadTestsFromTestCase(TestLoadTree)
])


if __name__ == '__main__':
    # run all defined TestCases for only this module if invoked directly
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(AllTestsTrees)

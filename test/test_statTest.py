import numpy as np
import pandas as pd

import os
import sys
sys.path.append('..')
from variationanalysis import score, statTest

import unittest

class TestStatTest(unittest.TestCase):
    
    def test_all_zeros(self):
        inp = np.array([[0]*np.random.randint(1e5)])
        self.assertEqual(statTest.Ttest(inp, save=False), (np.array([1]), np.array([0])))
        
    def test_symmetric_around_zero(self):
        inp = np.array([[-2.5,2.5]*np.random.randint(1e2)])
        self.assertEqual(statTest.Ttest(inp, save=False), (np.array([1]), np.array([0])))
        
    def test_all_pos(self):
        inp = np.array([[2]*np.random.randint(1e2)])
        self.assertEqual(statTest.Ttest(inp, save=False), (np.array([0]), np.array([np.inf])))
    
    def test_all_neg(self):
        inp = np.array([[-2]*np.random.randint(1e2)])
        self.assertEqual(statTest.Ttest(inp, save=False), (np.array([0]), np.array([-np.inf])))
    

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStatTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

    
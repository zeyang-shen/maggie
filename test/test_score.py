import numpy as np
import pandas as pd

import os
import sys
sys.path.append('..')
from variationanalysis import score, statTest

import unittest

class TestScoreDifference(unittest.TestCase):
    
    def test_all_zeros(self):
        inp = pd.DataFrame([[0]*10]*20)
        oup = score.compute_score_difference(inp)
        self.assertTrue((oup==0).all().all())
        
    def test_default_direction(self):
        inp = pd.DataFrame([[0]*10, [1]*10]*5)
        oup = score.compute_score_difference(inp)
        self.assertTrue((oup==-1).all().all())
        
    def test_reverse_direction(self):
        inp = pd.DataFrame([[0]*10, [1]*10]*5)
        oup = score.compute_score_difference(inp, reverse=True)
        self.assertTrue((oup==1).all().all())
    

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestScoreDifference)
    unittest.TextTestRunner(verbosity=2).run(suite)

    
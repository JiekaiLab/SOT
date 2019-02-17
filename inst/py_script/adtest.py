# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 10:33:27 2019

@author: lhlin
"""

from scipy.stats import anderson_ksamp
import numpy as np

def py_adtest(mat, lv):
    ds = set(lv)
    pv = []
    for i in np.arange(mat.shape[0]):
        pv.append(anderson_ksamp([mat[i,np.array(lv)==l] for l in ds]).significance_level)
    return pv

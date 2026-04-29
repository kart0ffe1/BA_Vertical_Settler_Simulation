# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 19:47:16 2026

@author: kaku03
"""

import numpy as np

def r_squared(exp, sim):
    enum = 0
    for i in range(0, len(exp)):
        enum += (exp[i]- sim[i])**2
    denom = 0
    for i in range(0, len(exp)):
        denom += (exp[i] - np.mean((exp)))**2
    return 1 - enum/denom

def nrmse(exp, sim, norm):
    enum = 0
    for i in range(0, len(exp)):
        enum += (exp[i] -  sim[i])**2
    return ((enum/len(exp))**0.5)/norm
    
    
        
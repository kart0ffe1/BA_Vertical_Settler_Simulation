# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 04:02:34 2026

@author: kaku03
"""

import numpy as np

class Simulation:
    def __init__(self, Exp):
        self.Exp = Exp
        self.delta_t = 0.5;
        self.time = None;
        self.h_c = None;
        self.h_s = None;
        self.h_p = None;
    
    def simulate(self):
        self.time = []
        self.h_c = []
        t = 0
        while t < self.Exp.t_i:
            self.time.append(t)
            self.h_c.append((self.Exp.H_0*1e-3 - self.Exp.V*1e-3*t + (self.Exp.V*1e-3*self.Exp.delta_h_i*1e-3)/self.Exp.psi_i*1e-3*(1 - np.exp(-self.Exp.psi_i*1e-3*t/(self.Exp.delta_h_i*1e-3))))*1e3)
            t += self.delta_t
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 23:36:00 2026

@author: kaku03
"""

import Experiment as exp
import Simulation as sim

file_dir = r'\\avt.rwth-aachen.de\home$\student\kaku03\.AVT-UserConfig\Documents\Messdaten\Hold-Up 0_3\457 mm Hold-Up 0_3.xlsx'
exp = exp.Experiment(file_dir)
exp.plot_experiment()

sim = sim.Simulation(exp)
sim.simulate()
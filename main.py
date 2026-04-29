# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 23:36:00 2026

@author: kaku03
"""

import Experiment as exp
import Simulation as sim
file_dir = r'\\avt.rwth-aachen.de\home$\student\kaku03\.AVT-UserConfig\Documents\Messdaten\Hold-Up 0_3\183 mm Hold-Up 0_3.xlsx'
exp = exp.Experiment(file_dir)
exp.plot_experiment()

sim1 = sim.Simulation(exp)
sim1.calc_v_0(starting_range=2)

#sim1.calc_phi_0_Pilhofer_Mewes()
sim1.calc_phi_0_Hartland_Kumar()

#sim1.calc_tau_0_exp()
sim1.calc_tau_0_phys()

sim1.simulate()
sim1.plot_Results()
print(sim1.err())

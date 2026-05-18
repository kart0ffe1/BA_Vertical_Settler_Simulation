# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 23:36:00 2026

@author: kaku03
"""

import Experiment as exp
import Simulation as sim
import Sed_Model as sed
import numpy as np

file_dir = r'\\avt.rwth-aachen.de\home$\student\kaku03\.AVT-UserConfig\Documents\Messdaten\Hold-Up 0_3\915 mm Hold-Up 0_3.xlsx'
exp = exp.Experiment(file_dir)
exp.plot_experiment()

sed = sed.Sed_Model(exp)
sed.calc_v_0(starting_range=3)
sed.calc_phi_0_Hartland_Kumar()
#sed.calc_phi_0_Pilhofer_Mewes()

sim1 = sim.Simulation(exp, sed)


#sim1.calc_tau_0_exp()
#sim1.calc_tau_0_phys()
sim1.calc_tau_0_num()
#sim1.set_parameters(9.2*1e-3, 0.847*1e-3, 24)    

sim1.simulate()
sim1.plot_Results()
print(sim1.err())

#h_c_t_f = np.interp(173.8, sim1.time, sim1.h_c)
#h_s_t_f = np.interp(173.8, sim1.time, sim1.h_s)

#abb = h_c_t_f - h_s_t_f
#print(abb)

#sim1.exp_Params(r'Export\183_val.xlsx')

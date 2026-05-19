# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 22:55:18 2026

@author: kaku03
"""
import pandas as pd
import matplotlib.pyplot as plt

class Batch_Experiment:
    def __init__(self, file_dir):
        self.df = pd.read_excel(file_dir, header=None)
        self.H_0 = self.df.iloc[1,0] * 1e-3
        self.Epsilon_0 = self.df.iloc[1,1]
        self.Epsilon_p = self.df.iloc[1,2]
        self.Mu_c = self.df.iloc[1,3] * 1e-3
        self.Mu_d  =self.df.iloc[1,4]*1e-3
        self.Rho_c = self.df.iloc[1,5]
        self.Delta_rho = self.df.iloc[1,6]
        self.sigma = self.df.iloc[1,7] *1e-3
        self.G = 9.81
        self.A_m = 2.15*1e-21
       
        hc_data = self.df.iloc[5:, [0,1]].dropna()
        hs_data = self.df.iloc[5:, [2,3]].dropna()
        self.H_c_time = hc_data.iloc[:,0].to_numpy() 
        self.H_c = hc_data.iloc[:,1].to_numpy()*1e-3
        self.H_s_time = hs_data.iloc[:,0].to_numpy()
        self.H_s = hs_data.iloc[:,1].to_numpy()*1e-3      
        self.T_f = max(self.H_s_time[-1], self.H_c_time[-1])
       
        
    def plot_experiment(self):        
        plt.figure(num='Experiment Data')
        plt.scatter(self.H_c_time, self.H_c, c='blue', label= 'Coalesence')
        plt.scatter(self.H_s_time, self.H_s, c='red', label= 'Sedimentation')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Height [m]')
        plt.show()
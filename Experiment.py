# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 22:55:18 2026

@author: kaku03
"""
import pandas as pd
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt

class Experiment:
    def __init__(self, file_dir):
        self.df = pd.read_excel(file_dir, header=None)
        self.df = self.df.dropna(how='all')
        self.H_0 = self.df.iloc[1,0]
        self.Epsilon_0 = self.df.iloc[1,1]
        self.Epsilon_p = self.df.iloc[1,2]
        self.Mu_c = self.df.iloc[1,3]
        self.Rho_c = self.df.iloc[1,4]
        self.Delta_roh = self.df.iloc[1,5]
        self.G = 9.81
       
        hc_data = self.df.iloc[5:, [0,1]].dropna()
        hs_data = self.df.iloc[5:, [2,3]].dropna()
        self.H_c_time = hc_data.iloc[:,0].to_numpy()
        self.H_c      = hc_data.iloc[:,1].to_numpy()
        self.H_s_time = hs_data.iloc[:,0].to_numpy()
        self.H_s      = hs_data.iloc[:,1].to_numpy()        
       
        self.v_0 = None
        self.phi_0 = None
        self.psi_i = None
        self.V = None
        self.t_i = None
        self.tau_0 = None
        self.delta_h_i = None
        
    def calc_v_0(self):
        self.v_0 = (self.H_s[1] - self.H_s[0])/(self.H_s_time[1] - self.H_s_time[0])
        
    def calc_phi_0(self):
        if self.v_0 is None:
            raise ValueError('v_0 has not been calculated')
        
        def f(x):
            return self.v_0*1e-3 - 12*self.Mu_c*1e-3/(0.53*self.Rho_c*x)*(-1 + np.sqrt(1 + (0.53*self.Rho_c*self.Delta_roh*self.G*x**3*(1 - self.Epsilon_0))/(108*(self.Mu_c*1e-3)**2*(1 + 4.56*self.Epsilon_0**0.73))))
        
        self.phi_0 = optimize.brentq(f, 1e-2, 1e-5)*1e3
       
    def calc_psi_i(self):
        max_incl = abs((self.H_c[1] - self.H_c[0])/(self.H_c_time[1]- self.H_c_time[0]))
        for i in range(1, len(self.H_c)-1):
            if abs((self.H_c[i+1] - self.H_c[i])/(self.H_c_time[i+1] -  self.H_c_time[i])) > max_incl:
                max_incl = abs((self.H_c[i+1] - self.H_c[i])/(self.H_c_time[i+1] -  self.H_c_time[i]))
                
        self.psi_i = max_incl
        
    def calc_t_i(self):
        if self.v_0 is None:
            raise ValueError('v_0 has not been calculated')
        if self.psi_i is None:
            raise ValueError('psi_i has not been calculated')
            
        def f(x):            
            if  self.psi_i*1e-3/(self.H_0*1e-3/x - self.v_0*1e-3/2 - (1 -self.Epsilon_p)*self.psi_i*1e-3/(2*self.Epsilon_p)) >= 1:
                return 1e10
            
            return (1 - self.Epsilon_0)/(1 - self.Epsilon_p)*self.H_0*1e-3 - (self.Epsilon_p*self.v_0*1e-3*x)/(2*(1 - self.Epsilon_p)) - (self.psi_i*1e-3*x)/2 - self.H_0*1e-3 + (self.H_0*1e-3/x - self.v_0*1e-3/2 - (1 -self.Epsilon_p)*self.psi_i*1e-3/(2*self.Epsilon_p))*x + self.psi_i*1e-3*x/(np.log(1 - self.psi_i*1e-3/(self.H_0*1e-3/x - self.v_0*1e-3/2 - (1 -self.Epsilon_p)*self.psi_i*1e-3/(2*self.Epsilon_p))))
        
        self.t_i = optimize.newton(f, self.H_c_time[-1]/2, maxiter=100)
        self.V = (self.H_0*1e-3/self.t_i - self.v_0*1e-3/2 - (1 -self.Epsilon_p)*self.psi_i*1e-3/(2*self.Epsilon_p))*1e3
        
    def calc_tau_0_exp(self):
        if self.v_0 is None: 
            raise ValueError('v_0 has not been calculated')
        if self.psi_i is None:
            raise ValueError('psi_i has not been calculated')
        if self.t_i is None:
            raise ValueError('t_i has not been calculated')
            
        self.tau_0 = (2*self.H_0*1e-3*self.Epsilon_p*(1 - self.Epsilon_0) - (self.v_0*1e-3*self.Epsilon_p + self.psi_i*1e-3*(1 - self.Epsilon_p))*self.t_i)/(3*self.psi_i*1e-3*(1 - self.Epsilon_p))
   
    def calc_delta_h_i(self):
        if self.v_0 is None:
            raise ValueError('v_0 has not been calculated')
        if self.psi_i is None:
            raise ValueError('psi_i has not been calculated')
        if self.t_i is None:
            raise ValueError('t_i has not been calculated')
        
        self.delta_h_i = ((1 - self.Epsilon_0)/(1 - self.Epsilon_p)*self.H_0*1e-3 - (self.v_0*1e-3*self.t_i)/(2*(1 - self.Epsilon_p)) - (self.psi_i*1e-3*self.t_i)/(2*self.Epsilon_p))*1e3
        
    def calc_params(self):
        self.calc_v_0()
        self.calc_phi_0()
        self.calc_psi_i()
        self.calc_t_i()
        self.calc_tau_0_exp()
        self.calc_delta_h_i()
        
    def plot_experiment(self):
        self.calc_params()
        
        plt.figure(num='Experiment Data')
        plt.scatter(self.H_c_time, self.H_c, c='blue', label= 'Coalesence')
        plt.scatter(self.H_s_time, self.H_s, c='red', label= 'Sedimentation')
        plt.vlines(self.t_i, ymin= self.H_s[0], ymax=self.H_c[0], linestyles='dashed', colors='black', label= 't_i')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Height [mm]')
        plt.show()
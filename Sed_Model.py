# -*- coding: utf-8 -*-
"""
Created on Mon May 18 10:56:06 2026

@author: kaku03
"""

from scipy import optimize
import helper_functions as hlp
from matplotlib import pyplot as plt
import numpy as np

class Sed_Model:
    def __init__(self, Exp):
        self.Exp = Exp
        self.v_0 = None
        self.phi_0 = None
    
    def calc_v_0(self, starting_range = 3, plot= True):
        reg_val= []
        func_val = []   
         
        for i in range(0,starting_range):
            reg_val.append(self.Exp.H_s[i])
        
        def f(t, v_0, b):
            return t*v_0 + b
        
        popt, covariance = optimize.curve_fit(f, self.Exp.H_s_time[0:len(reg_val)], reg_val)        
        v_0, b = popt
        
        for i in range (0, len(reg_val)):
            func_val.append(f(self.Exp.H_s_time[i], v_0, b)) 
        
        r_sqr = hlp.r_squared(reg_val, func_val)
        
        for i in range(starting_range, len(self.Exp.H_s)):
            reg_val_new = reg_val.copy()
            func_val_new = []
            
            reg_val_new.append(self.Exp.H_s[i])
            popt, covariance = optimize.curve_fit(f, self.Exp.H_s_time[0:len(reg_val_new)], reg_val_new)
            v_0_new, b_new = popt
            
            for j in range (0, len(reg_val_new)):
                func_val_new.append(f(self.Exp.H_s_time[j], v_0_new, b_new))
                
            r_sqr_new = hlp.r_squared(reg_val_new, func_val_new)
            
            if abs(v_0) > abs(v_0_new) and r_sqr > r_sqr_new:
                break
           
            reg_val = reg_val_new.copy()
            func_val = func_val_new.copy()
            v_0 = v_0_new
            b  = b_new
            
        self.v_0 = v_0
        
        if plot:
            plt.figure(num='Linear Regression Sedimentation')
            plt.scatter(self.Exp.H_s_time, self.Exp.H_s, c='red', label='Sedimentation Measured')
            plt.plot(self.Exp.H_s_time[0:len(func_val)],func_val, c='blue', linestyle='dashed', label='Sedimentation Regression')
            plt.vlines(self.Exp.H_s_time[0], 0, self.Exp.H_s[-1], colors='black', linestyles='dashdot', label='Start of Regression')
            plt.vlines(self.Exp.H_s_time[len(func_val) - 1], 0, self.Exp.H_s[-1], colors='black', linestyles='solid', label='End of Regression')
            plt.xlabel('Time [s]')
            plt.ylabel('Height [m]')
            plt.legend()
            plt.show()
        
    def calc_phi_0_Hartland_Kumar(self):
        if self.v_0 is None:
            raise ValueError('v_0 has not been calculated')

        def phi_0(phi_0):
            return self.v_0 - 12*self.Exp.Mu_c/(0.53*self.Exp.Rho_c*phi_0)*(-1 + np.sqrt(1 + (0.53*self.Exp.Rho_c*self.Exp.Delta_rho*self.Exp.G*phi_0**3*(1 - self.Exp.Epsilon_0))/(108*(self.Exp.Mu_c)**2*(1 + 4.56*self.Exp.Epsilon_0**0.73))))
            
        self.phi_0 = optimize.brentq(phi_0, 1e-5, 1e-2)
        
    def calc_phi_0_Pilhofer_Mewes(self):
        if self.v_0 is None:
            raise ValueError('v_0 has not been calculated')
        
        
        def Ar(phi_0):
            Ar = (
                self.Exp.G*phi_0**3*self.Exp.Delta_rho*self.Exp.Rho_c
                /(self.Exp.Mu_c**2)
            )            
            return Ar
        
        def K_HR():
            K_HR = (
                3*(self.Exp.Mu_c + self.Exp.Mu_d)
                /(2*self.Exp.Mu_c + 3*self.Exp.Mu_d)
            )            
            return K_HR
        
        def zeta():
            zeta = (
                5*K_HR()**(-3/2)*(self.Exp.Epsilon_0
                /(1 - self.Exp.Epsilon_0))**0.45
            )            
            return zeta
        
        def q():
            q = (
                (1 - self.Exp.Epsilon_0)/(2*self.Exp.Epsilon_0
                *K_HR())*np.exp(2.5*self.Exp.Epsilon_0
                /(1 - (0.61*self.Exp.Epsilon_0)))
            )
            return q
        
        def Re_inf(phi_0):
            Re_inf = (
                9.72*((1 + 0.01*Ar(phi_0))**(4/7) - 1)    
            ) 
            return Re_inf
        
        def cw(phi_0):
            cw = (
                Ar(phi_0)/(6*Re_inf(phi_0)**2) 
                - 3/(K_HR()*Re_inf(phi_0))    
            )
            return cw
        
        def v_rs():
            return self.v_0/(1 - self.Exp.Epsilon_0)
        
        def Re_s(phi_0):
            return self.Exp.Rho_c*v_rs()*phi_0/self.Exp.Mu_c  
        
        def v_Stokes(phi_0):
            v_Stokes = (
                self.Exp.G*phi_0**2*self.Exp.Delta_rho
                /(18*self.Exp.Mu_c)
            )
            return v_Stokes
        
        def Re_s_iteration(phi_0):
            a = (
                3*q()*self.Exp.Epsilon_p/(cw(phi_0)*zeta()
                *(1 - self.Exp.Epsilon_0))    
            )
            
            b = (
                ((1 + Ar(phi_0)*(cw(phi_0)*zeta()*
                (1 - self.Exp.Epsilon_0)**3)/(54*q()**2
                *self.Exp.Epsilon_0**2))**(1/2) - 1)
            )
            return a*b - Re_s(phi_0)
        
        res = optimize.brentq(Re_s_iteration, 1e-8, 10e-3, full_output= True, xtol=2e-12)
        self.phi_0 = res[0]
    
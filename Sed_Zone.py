# -*- coding: utf-8 -*-
"""
Created on Mon May 18 11:36:19 2026

@author: kaku03
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize as optimize

class Sed_Zone:
    def __init__(self):
        self.Roh_c = None
        self.Delta_roh = None
        self.Mu_c = None
        self.G = 9.81
        
        self.delta_x = 0.01
        self.q_d = None
        self.q_e =  None
        #self.R = self.q_d/self.q_e
        self.R = None
        self.epsilon_0 = None
        self.epsilon_s = None
        self.tau_0 = None
        self.phi_0 = None
        
        self.x = None
        self.epsilon = None
        self.tau = None
        self.phi = None
        self.H_s = None
        
    def calc_tau(self,tau_0):
        tau = tau_0*self.phi[-1]/self.phi_0 
        return tau
    
    def calc_phi(self, turb=False):
        if not turb:
            phi_sqr = (
                18*self.Mu_c*self.q_d/(self.Delta_roh*self.G)*(1 + (self.R + 3.56)
                *self.epsilon[-1] + 4.56*(self.R - 1)*self.epsilon[-1]**2)
                /(self.epsilon[-1]*(1 - self.epsilon[-1])**2)
            )
            
            return phi_sqr**(1/2)
        
        else:
            phi = (
                1.59*self.Roh_c*self.q_d**2/(4*self.Delta_roh*self.G)
                *(1 + 4.56*self.epsilon[-1])*(1 + (self.R - 1)*self.epsilon[-1])**2
                /(self.epsilon[-1]**2*(1 - self.epsilon[-1])**3)
            )
            
            return phi
    
    def calc_epsilon(self, turb=False):
        if not turb:
            def f(x):
                f = (
                    self.x[-1] - 3*self.q_d*self.tau[-1]*np.log(
                    (x/self.epsilon_0)**(self.R + 5.56)
                    *((1 - self.epsilon_0)/(1 - x))**2
                    *((1 + 4.56*self.epsilon_0)/(1 + 4.56*x))**4.56
                    *((1 + (self.R - 1)*self.epsilon_0)
                    /(1 + (self.R - 1)*x))**(self.R - 1))
                    + 3*self.q_d*self.tau[-1]*(1/self.epsilon_0 - 1/x)
                    )
        
                return f
            
            # val = np.linspace(0.01, 0.99, num = 100)
            # test = []
            
            # for x in val:
            #     test.append(f(x))
                
            # plt.plot(val, test)
            # plt.show()
            
            epsilon = optimize.brentq(f, self.epsilon_0, self.epsilon_s)
           
            
        else:
            def f(x):
                f = (
                    self.x[-1] - 6*self.q_d*self.tau[-1]*np.log(
                    (x/self.epsilon_0)**(2*self.R + 5.56)
                    *((1 - self.epsilon_0)/(1 - x))**3
                    *((1 + 4.56*self.epsilon_0)/(1 + 4.56*x))**4.56
                    *((1 + (self.R - 1)*self.epsilon_0)
                    /(1 + (self.R - 1)*x))**(2*self.R - 2))
                    + 12*self.q_d*self.tau[-1]*(1/self.epsilon_0 - 1/x)
                )
                
                return f
            
            epsilon = optimize.brentq(f, self.epsilon_0, self.epsilon_s)
        return epsilon        
            
        
    
    def simulate(self, turb=False):
        self.x = []
        self.x.append(0)
        
        self.phi = []
        self.phi.append(self.phi_0)
        
        self.tau = []
        self.tau.append(self.tau_0)
        
        self.epsilon = []
        self.epsilon.append(self.calc_epsilon(turb))        
        
        while self.epsilon[-1] < self.epsilon_s:
            self.x.append(self.x[-1] + self.delta_x)
            self.phi.append(self.calc_phi(turb))
            self.tau.append(self.calc_tau(self.tau_0))
            try:
                self.epsilon.append(self.calc_epsilon(turb))
            except:
                self.epsilon.append(self.epsilon_s)
                
                # print(self.x[-1])
                # print(self.epsilon[-1])
                # def f(x):
                #     f = (
                #         self.x[-1] - 6*self.q_d*self.tau[-1]*np.log(
                #         (x/self.epsilon_0)**(2*self.R + 5.56)
                #         *((1 - self.epsilon_0)/(1 - x))**3
                #         *((1 + 4.56*self.epsilon_0)/(1 + 4.56*x))**4.56
                #         *((1 + (self.R - 1)*self.epsilon_0)
                #         /(1 + (self.R - 1)*x))**(2*self.R - 2))
                #         + 12*self.q_d*self.tau[-1]*(1/self.epsilon_0 - 1/x)
                #     )
                    
                #     return f
                                
                # val = np.linspace(self.epsilon_0, self.epsilon_s, num = 100)
                # test = []
                # for x in val:
                #     test.append(f(x))
                
                # plt.plot(val, test)
                # plt.show()
                # plt.plot(self.x[0:-1], self.epsilon)
                # plt.show()
                # raise ValueError('BrentQ failed')
                
            
            
        self.H_s = np.interp(self.epsilon_s, self.epsilon, self.x)
    
    def plot_Results(self, dim = True):
        if dim:        
            plt.figure(num='Simulation')
            plt.plot(self.x, self.epsilon, c='blue', label='Hold Up')
            plt.vlines(self.H_s, 0, 1, linestyles='dashed', colors='black', label='H_s')
            plt.legend()
            plt.xlabel('h [m]')
            plt.ylabel('Hold Up [-]')
            plt.show()
        
        else:
            plt.figure(num='Dimensionless Simulation')
            
            x_dimless = []
            epsilon_dimless = []
            
            for i in range(0, len(self.x)):
                x_dimless.append(self.x[i]/self.H_s)
                epsilon_dimless.append(self.epsilon[i]/self.epsilon_s)
            
            plt.plot(x_dimless, epsilon_dimless, c='blue', label='Hold Up')
            plt.legend()
            plt.xlabel('h/H_s [-]')
            plt.ylabel('ɛ/ɛ_s [-]')
            plt.show()
        
        
        
        
            
        
            
        
    
    
    
        
            
            

       
        
        
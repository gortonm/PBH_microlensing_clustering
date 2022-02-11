#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 17:49:13 2022

@author: ppxmg2
"""
import numpy as np
import halomodel_optimised as hm
import os

m_pbh = 1
n_cl = 10**5
setup = hm.Halomodel(m_pbh=m_pbh, n_cl=n_cl)

speed_conversion = 1.022704735e-6
c, G = 2.99792458e5 * speed_conversion, 4.30091e-3 * speed_conversion**2 

def trapezium(x, y):
    area = 0
    for i in range(len(x) - 1):
        area += (x[i+1] - x[i]) * (y[i] + 0.5*(x[i+1]-x[i])*(y[i+1]-y[i]))
    return area
    
def dGamma_integrand(setup, d_L, t_hat):
    integrand = 512 * setup.rho_0 * (setup.r_0**2 + setup.a_0**2) * G**2 * m_pbh * np.exp(-(2 * setup.einstein_radius(d_L) / (t_hat * setup.v_c))**2) * d_L**2 * (setup.d_s - d_L)**2 * (d_L**2 + (setup.a*d_L) + (setup.r_0**2 + setup.a_0**2))**(-1) / (setup.d_s**2 * t_hat**4 * setup.v_c**2 * c**4)
    return integrand

dgammas_sum_10000 = []
path_dir = f'{os.getcwd()}'

t_hat_values = np.arange(1, 5000, 1.) / 365.25
for m_pbh in [1, 10, 100]:
    for step in [0.5]:
        setup = hm.Halomodel(m_pbh=m_pbh, n_cl=n_cl)
        d_Ls = np.arange(0., setup.d_s, step)
        dgammas = []          
        
        for t_hat in t_hat_values:
            dgammas.append(trapezium(d_Ls, dGamma_integrand(setup, d_Ls, t_hat)))
        
        np.savetxt(path_dir + '/data_files/event_duration_dist_step={:.2f}_mpbh={:.2f}'.format(step, np.log10(m_pbh)) + '.txt', (np.array(t_hat_values) * 365.25, np.array(dgammas)), delimiter = ', ')
    
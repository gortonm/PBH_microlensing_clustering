#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 15:23:27 2022
@author: ppxmg2
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from expected_events_discrete_clustered import n_ex
import halomodel_optimised as hm

new_RE = False

if new_RE:
    append = 'newRE'
    
else:
    append = 'oldRE'


# Range of PBH masses to consider
#m_pbhs = 10**np.arange(0., 1.1, 1)
m_pbhs = 10**np.array([0.])

# Range of numbers of PBHs per cluster to consider
#n_cls = 10**np.arange(7., 2.9, -1.)
n_cls = 10**np.array([6.])

# total number of realisations
n_realisations = 1000


"""Calculate the mean number of expected events"""
for m_pbh in m_pbhs:
   
    for j, n_cl in enumerate(n_cls):
        
        setup = hm.Halomodel(n_cl=n_cl, m_pbh=m_pbh)
        dL_lim = setup.dL_lim
        
        x_min = []
        n_ex_mean = []
        n_ex_poisson = []
        
        x_min_min = 0.05
        i_min = 0
        
        n_realisations_small_x = 0
        n_samples_small_x = 0
        n_samples = 0
        dLs_length = 0
       
        
        for i in range(n_realisations):
            
            np.random.seed = int(m_pbh*n_cl)
            
            filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + str(i)
            d_L = np.loadtxt(filepath + append + '_dL.txt', delimiter=',')
            v = np.loadtxt(filepath + append + '_v.txt', delimiter=',')

            # calculate number of expected events
            n_ex_mean.append(n_ex(d_L, v, m_pbh, n_cl, eff=True, poisson=False))
            n_ex_poisson.append(n_ex(d_L, v, m_pbh, n_cl, eff=False, poisson=True))
            
            x_min.append(min(d_L) / 50e3)
            
            if min(d_L) / 50e3 < x_min_min:
                x_min_min = min(d_L) / 50e3
                i_min = i
            
            
            
            
            if min(d_L) / 50e3 < 0.01:
                print('i = {:.0f}'.format(i) + ', x={:.5f}'.format(min(d_L)/50e3))
            
            
            n_samples += len(d_L)
            n_samples_small_x += len(d_L[d_L < dL_lim])
            if min(d_L) < dL_lim:
                n_realisations_small_x += 1
                
        
        print('N_cl=1e', np.log10(n_cl))
        print('M_PBH=1e', np.log10(m_pbh))

        print('i_min = ', i_min)
        print(n_realisations_small_x)

        print(n_samples_small_x / n_samples)

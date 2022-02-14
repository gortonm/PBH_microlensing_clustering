#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 15:08:32 2022

@author: ppxmg2
"""

import numpy as np
import os
from tE_distribution_producevalues_optimised_v2 import produce_values
from expected_events_discrete_clustered_optimised import n_ex

np.random.seed = 8022022

# Range of PBH masses to consider
#m_pbhs = 10**np.arange(-1, 3.1, 1)
m_pbhs = 10**np.arange(0, 2.1, 1)

# Range of numbers of PBHs per cluster to consider
#n_cls = 10**np.arange(3, 7.1, 1.)
n_cls = 10**np.arange(3, 5.1, 1.)

# Number of realisations for each PBH mass and cluster size
n_realisations = 20

# Set parameters for the standard halo model with the LMC at a distance of 
# 50 kpc and a circular velocity of the Sun of 220 km/s
d_s, v_c = 50, 220

"""Produce a sample of cluster line of sight distances and speeds"""
for m_pbh in m_pbhs:
    for n_cl in n_cls:
                
        for i in range(n_realisations):
            d_L, v = produce_values(n_cl, m_pbh, d_s, v_c)
            
            # Save cluster distances and speeds
            filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + str(i)
            np.savetxt(filepath + '_dL.txt', d_L)
            np.savetxt(filepath + '_v.txt', v)
            
            if i % 5 == 0:
                print(i)
            
"""Calculate the number of expected events, given a sample of cluster distances and speeds"""
for m_pbh in m_pbhs:
    
    for n_cl in n_cls:
        n_ex_EROS_efficiency = 0
        n_ex_EROS_efficiency_blendingcorrection = 0
        n_ex_perfect_efficiency = 0
        
        for i in range(n_realisations):
            filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + str(i)
            d_L = np.loadtxt(filepath + '_dL.txt', delimiter=',')
            v = np.loadtxt(filepath + '_v.txt', delimiter=',')
            
            # calculate number of expected events
            n_ex_EROS_efficiency += n_ex(d_L, v, m_pbh, n_cl, eff=True)    # EROS-2 efficiency curve
            n_ex_EROS_efficiency_blendingcorrection += n_ex(d_L, v, m_pbh, n_cl, eff=True, blendingcorrection = True)    # EROS-2 efficiency curve, with average blending correction
            n_ex_perfect_efficiency += n_ex(d_L, v, m_pbh, n_cl, eff=False)    # perfect efficiency
            
        np.savetxt(filepath + '_n_ex_EROS.txt', n_ex_EROS_efficiency)
        np.savetxt(filepath + '_n_ex_EROS_blendingcorrection.txt', n_ex_EROS_efficiency_blendingcorrection)
        np.savetxt(filepath + '_n_ex_perfect.txt', n_ex_perfect_efficiency)

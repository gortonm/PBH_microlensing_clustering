#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 15:23:27 2022

@author: ppxmg2
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.special import factorial

np.random.seed = 8022022

# return probability mass function for a Poisson distribution
def poisson_pmf(k, lam):
    return np.power(lam, k) * np.exp(-lam) / factorial(k)


# Range of PBH masses to consider
m_pbhs = 10**np.arange(1., 1.1, 1)

# Range of numbers of PBHs per cluster to consider
n_cls = 10**np.arange(3, 7.1, 1.)

# Number of realisations for each PBH mass and cluster size
n_realisations = 1000


for m_pbh in m_pbhs:
    
    # mean value of the number of expected events, from smoothly-distributed PBHs of mass m_pbh
    """TEMPORARY VALUE - NEEDS CALCULATING"""
    mean_smooth = 25.4

    
    colors = ['k', 'r', 'orange']
    bin_spacing = 1.
    bins = np.arange(0, 1000, bin_spacing)
    plt.figure(figsize=(10, 7))
    
    i = 0
    for n_cl in n_cls:
        color = colors[i]
        
        # load file for number of events in a particular realisation
        filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh))
        n_ex_EROS_efficiency = np.loadtxt(filepath + '_n_ex_EROS.txt')
        n_ex_EROS_efficiency_blendingcorrection = np.loadtxt(filepath + '_n_ex_EROS_blendingcorrection.txt')
        n_ex_perfect_efficiency = np.loadtxt(filepath + '_n_ex_perfect.txt')
      
        plt.hist(n_ex_EROS_efficiency, bins=bins, density=True, color=colors[i], histtype='step', label='$N_{cl} = $'+'$10^{:.0f}$'.format(np.log10(n_cl)))
        plt.hist(n_ex_EROS_efficiency_blendingcorrection, bins=bins, density=True, color=colors[i], histtype='step', label='With average EROS-2 blending correction: $N_{cl} = $'+'$10^{:.0f}$'.format(np.log10(n_cl)))        
        plt.hist(n_ex_perfect_efficiency, bins=bins, density=True, color=colors[i], histtype='step', linestyle='dotted', label='$\epsilon(t_E) = 1$, $N_{cl} = $'+'$10^{:.0f}$'.format(np.log10(n_cl)))
     
        plt.xlim(0, 100)
        
        i += 1
    
    # add Poisson distribution plot for mean value from smooth case
    plt.plot(poisson_pmf(np.arange(0, 101, 1.), lam=mean_smooth, color='k', label='Smooth PBH distribution'))
    plt.legend()
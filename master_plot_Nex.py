#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 15:23:27 2022

@author: ppxmg2
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from scipy.special import factorial


#Specify the plot style
mpl.rcParams.update({'font.size': 16,'font.family':'serif'})
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['xtick.top'] = False
mpl.rcParams['ytick.right'] = False
mpl.rcParams['font.family'] = 'serif'
mpl.rc('text', usetex=True)

mpl.rcParams['legend.edgecolor'] = 'lightgrey'

np.random.seed = 8022022

# return probability mass function for a Poisson distribution
def poisson_pmf(k, lam):
    return np.power(lam, k) * np.exp(-lam) / factorial(k)


# Range of PBH masses to consider
m_pbhs = 10**np.arange(0., 0.1, 1)

# Range of numbers of PBHs per cluster to consider
n_cls = 10**np.arange(6, 8.1, 1.)

# Number of realisations for each PBH mass and cluster size
n_realisations = 1000

# Approximate mean value for smooth case
n_cl, m_pbh = 10**3, 1.
filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh))
mean_smooth = 24.6
mean_smooth_perfect = 60.1

for m_pbh in m_pbhs:
    
    # mean value of the number of expected events, from smoothly-distributed PBHs of mass m_pbh
    """TEMPORARY VALUE - NEEDS CALCULATING"""

    colors = ['darkorange', 'r', 'saddlebrown', 'yellow']
    bin_spacing = 1
    bins = np.arange(0, 1000, bin_spacing)
    plt.figure(figsize=(5, 4))
    
    # add Poisson distribution plot for mean value from smooth case
    plt.plot(poisson_pmf(np.arange(0, 101, 1.), lam=mean_smooth), color='k', linestyle='dashed')
    plt.plot(poisson_pmf(np.arange(0, 101, 1.), lam=mean_smooth_perfect), color='k', linestyle='dotted')
    
    i = 0
    for n_cl in n_cls:
        color = colors[i]
        
        # load file for number of events in a particular realisation
        filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh))
        n_ex_EROS_efficiency = np.loadtxt(filepath + 'n_ex_EROS.txt')
        
        print(len(n_ex_EROS_efficiency))
        
        #n_ex_EROS_efficiency_blendingcorrection = np.loadtxt(filepath + '_n_ex_EROS_blendingcorrection.txt')
        n_ex_perfect_efficiency = np.loadtxt(filepath + 'n_ex_perfect.txt')
      
        plt.hist(n_ex_EROS_efficiency, bins=bins, density=True, color=colors[i], histtype='step', label='$N_\mathrm{cl} = $'+'$10^{:.0f}$'.format(np.log10(n_cl)))
        #plt.hist(n_ex_EROS_efficiency_blendingcorrection, bins=bins, density=True, color=colors[i], histtype='step', label='With average EROS-2 blending correction: $N_{cl} = $'+'$10^{:.0f}$'.format(np.log10(n_cl)))        
        #plt.hist(n_ex_perfect_efficiency, bins=bins, density=True, color=colors[i], histtype='step', linestyle='dotted', label='$\epsilon(t_E) = 1$, $N_{cl} = $'+'$10^{:.0f}$'.format(np.log10(n_cl)))
        plt.hist(n_ex_perfect_efficiency, bins=bins, density=True, color=colors[i], histtype='step', linestyle='dotted')
                 
        i += 1
        
    
    plt.xlim(0, 150)
    plt.xlabel('$N_\mathrm{obs}$')
    plt.ylabel('$P(N_\mathrm{obs})$')
    plt.tight_layout()
    plt.legend()

    # save figure without legend
    plt.savefig('N_exp_mpbh={:.2f}'.format((m_pbh)) + 'Msun_binspacing={:.0f}'.format(bin_spacing) + '.pdf')
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 14:03:30 2022

@author: ppxmg2
"""

import numpy as np
import os

# Boolean controlling whether to use EROS-2 efficiency function and exposure
EROS2_eff = False

# Boolean controlling whether to set PBH cluster radius to 10 pc
set_rcl_10 = True

append = ""
if set_rcl_10:
    append += "_rcl10_"
if EROS2_eff:
    append += "_EROS2_"


# Set number of PBHs per cluster and PBH mass (in solar masses)
n_cl, m_pbh = 1e3, 1e2
n_realisations = 300000

# Find probability corresponding to a given value of N_obs
def prob_n_obs(n_obs_data, n_obs):
    
    # Create histogram
    bins = np.arange(-0.5, max(n_obs_data)+0.51, 1)
    #bins = np.arange(0., max(n_obs_data)+0.51, 1)

    probability, bin_edges = np.histogram(n_obs_data, bins=bins, density=True)

    # Since bins are centred at integer values, and the first bin is centred at zero,
    # the index of the bin matches that of the integer value it corresponds to
    if n_obs < len(probability):
        return probability[n_obs]
    else:
        return 0

def f_pbh_proability_function(f_pbh):
    f_pbh = round(f_pbh, 4)
    
    n_obs_data = np.loadtxt(f"{os.getcwd()}" + "/simulated_data_constraints/N_cl/{0:.2f}".format(np.log10(n_cl)) + "/M_PBH/{0:.2f}/".format(np.log10(m_pbh)) + 'n_ex_corrected_fpbh={0:.4f}'.format(f_pbh) + append + '_nsamp={0:.1e}'.format(n_realisations) + '.txt')

    
    prob_new = prob_n_obs(n_obs_data, 0)
    print('f_PBH = ', f_pbh)
    print('P(N) = ', prob_new)
    print(prob_new - 0.05)
    return prob_new - 0.05
    
#f_pbhs=np.arange(0.020, 0.03, 0.001)
f_pbhs = np.arange(0.0200, 0.022, 0.0001)

#f_pbhs = np.arange(0.09, 0.101, 0.001)
#f_pbhs = np.arange(0.095, 0.09, 0.0001)
#f_pbhs = np.arange(0.094, 0.0976, 0.0001)

for f_pbh in f_pbhs:
    
    n_obs_data = np.loadtxt(f"{os.getcwd()}" + "/simulated_data_constraints/N_cl/{0:.2f}".format(np.log10(n_cl)) + "/M_PBH/{0:.2f}/".format(np.log10(m_pbh)) + 'n_ex_corrected_fpbh={0:.4f}'.format(f_pbh) + append + '_nsamp={0:.1e}'.format(n_realisations) + '.txt')

    print('f_PBH = {:.4f}, P(N=0) = {:.4f}'.format(f_pbh, prob_n_obs(n_obs_data, 0)))
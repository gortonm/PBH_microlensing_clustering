#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:12:28 2022

@author: ppxmg2
"""
from generate_results_updated_pcl_v2 import calculate_mass_cone_analytic
import numpy as np
import os

# Boolean controlling whether to use EROS-2 efficiency function and exposure
EROS2_eff = True

# Boolean controlling whether to set PBH cluster radius to 10 pc
set_rcl_10 = True

n_realisations = 10000

d_s = 50e3
x_max = 0.1

def convert_to_array(x):
    """
    Convert a scalar Float to a single-valued array, or do nothing if passed an array
    Parameters
    ----------
    x : Numpy array of type Float, or Float
    Returns
    -------
    Numpy array of type Float
        If passed scalar of type Float, return that value as a single-valued Numpy array of type Float.
        Otherwise, return the input array.
    """
    return np.array([x]) if isinstance(x, np.float64) else x



append = ""
if set_rcl_10:
    append += "_rcl10_"
if EROS2_eff:
    append += "_EROS2_"

if EROS2_eff:
    m_pbh = 1
    n_cls = 10**np.arange(6, 8.1, 1)
    
else:
    m_pbh = 1000
    n_cls = 10**np.arange(3, 5.1, 1)

for n_cl in n_cls:
    print('N_cl = {:.0e}'.format(n_cl))
    print('M_PBH / M_\odot = {:.0f}'.format(m_pbh))
    n_clusters_mean_extended = 0
    n_clusters_analytic = calculate_mass_cone_analytic() / (m_pbh * n_cl)
    # Load data from csv file
    
    for r in range(0, n_realisations):
        filename = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + "cluster_data_update_pcl_" + str(r) + append + ".csv"    
        
        data = np.genfromtxt(filename, delimiter=',')
        
        if len(data) == 0:
            n_clusters_mean_extended += 0
        else:
            dL = np.transpose(data)[0]
            dL_truncated = dL[dL < x_max * d_s]
            n_clusters_mean_extended += len(dL_truncated) / n_realisations
    
    print('Number of clusters in extended sample / Number of clusters in microlensing cone = ', n_clusters_mean_extended / n_clusters_analytic)


print(calculate_mass_cone_analytic(x_max))
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
import csv


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

mpl.rcParams['legend.edgecolor'] = 'inherit'


def find_quantiles(x, minval=0.16, maxval=0.84):
    x = np.sort(x)
    return x[int(len(x) * minval)-1], x[int(len(x) * maxval)-1]


def error_bars(x):
    # calculate upper and lower 68% and 95% percentiles
    x_68_lower, x_68_upper = find_quantiles(x) 
    x_95_lower, x_95_upper = find_quantiles(x, minval=0.025, maxval=0.975)
    
    # calculate arrays of 68% and 95% upper and lower error bars
    return [[np.mean(x) - x_68_lower], [x_68_upper - np.mean(x)]], [[np.mean(x) - x_95_lower], [x_95_upper - np.mean(x)]]


def load_constraints(Tisserand=True):
    """
    Load constraints on the PBH fraction from EROS-2, from Fig. 11 of Tisserand et al. (2007) (green curve)

    Returns
    -------
    m_pbh : List of type Float
        PBH mass (assumed monochromatic).
    f_pbh : List of type Float
        Constraint on the dark matter fraction in PBHs at a given (monochromatic)
        PBH mass m_pbh.

    """
    if Tisserand:
        name = 'tisserand'
    else:
        name = 'green_reproduction'
    file_constraints_CSV = open('./data_files/eros2_constraints_' + name + '.csv')
    constraints_CSV = csv.reader(file_constraints_CSV)
    list_constraints_CSV = list(constraints_CSV)
        
    m_pbh = []
    f_pbh = []
    for col in list_constraints_CSV[1:]:
        m_pbh.append(float(col[0]))
        f_pbh.append(float(col[1]))
        
    return m_pbh, f_pbh



# Range of PBH masses to consider
m_pbhs = 10**np.arange(0., 1.51, 0.5)

# Range of numbers of PBHs per cluster to consider
n_cls = 10**np.arange(3, 8.1, 1.)

# Number of realisations for each PBH mass and cluster size
n_realisations = 1000

# colours to plot for different cluster sizes
colours = ['darkturquoise', 'b', 'm', 'darkorange', 'r', 'saddlebrown']



plt.figure(figsize=(8,4))

""" Plot smooth constraint from EROS-2 (Fig. 15 of Tisserand et al. 2007 """
m_pbh_smooth, f_pbh_smooth = load_constraints()
plt.plot(m_pbh_smooth, f_pbh_smooth, color='lime')

""" Plot smooth constraint from EROS-2 (Fig. 1 of Green 2016 """
m_pbh_smooth, f_pbh_smooth = load_constraints(Tisserand=False)
plt.plot(m_pbh_smooth, f_pbh_smooth, linestyle='dotted', color='k')



for j, n_cl in enumerate(n_cls):
    
    # displace positions of different cluster sizes on the x-axis
    dx = 0.035 * (-len(n_cls)/2 + j*1.2)
    
    legend_list = []
    
    for m_pbh in m_pbhs:
        
        # load file for number of events in a particular realisation
        filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh))
        n_ex_EROS_efficiency = np.loadtxt(filepath + 'n_ex_EROS.txt')
        
        
        f_pbhs = np.ones(len(n_ex_EROS_efficiency))
        for i in range(len(n_ex_EROS_efficiency)):
            if n_ex_EROS_efficiency[i] > 3:
                f_pbhs[i] = (3 / n_ex_EROS_efficiency[i])
        
        # calculate arrays of 68% and 95% upper and lower error bars
        err_68, err_95 = error_bars(f_pbhs)
                
        # plot error bars and mean value
        plt.errorbar(np.log10(m_pbh) + dx, np.mean(f_pbhs), yerr = err_68, marker='x', color=colours[j], elinewidth=1, capsize=2)
        plt.errorbar(np.log10(m_pbh) + dx, np.mean(f_pbhs), yerr = err_95, marker='x', color=colours[j], elinewidth=0.5, capsize=1)
    
        print(n_cl)
        print(m_pbh)
        print(err_68)
        
    point = plt.plot(np.log10(m_pbh) + dx, np.mean(f_pbhs), color = colours[j], marker='x', linestyle='none', label='$10^{:.0f}$'.format(np.log10(n_cl)))

plt.legend(title='$N_\mathrm{cl}$')



plt.xlim(-0.4, 2.1)
plt.ylim(0., 1.01)

#plt.yscale('log')
plt.xlabel('$\log_{10}(M_\mathrm{PBH}/M_\odot)$')
plt.ylabel('$f_\mathrm{PBH}$')
plt.tight_layout()
plt.savefig('f_PBH_MPBH.pdf')

            
            
    
""" Test quantile finding method """
sample = np.random.uniform(size=1000000)
print(find_quantiles(sample))
print(find_quantiles(sample, minval=0.025, maxval=0.975))


err_68, err_95 = error_bars(sample)
plt.errorbar(0., np.mean(sample), yerr = err_68, marker='x', color='r', elinewidth=1, capsize=2)
plt.errorbar(0., np.mean(sample), yerr = err_95, marker='x', color='r', elinewidth=0.5, capsize=1)



""" Test matplotlib errorbar method """
lower_bound_1 = np.array([4, 7])
lower_bound_2 = np.array([2, 10])
mean = 5

err_68 = [[mean - lower_bound_1[0]], [lower_bound_1[1] - mean]]
err_95 = [[mean - lower_bound_2[0]], [lower_bound_2[1] - mean]]

plt.errorbar(0., mean, yerr = err_68, marker='x', color='r', elinewidth=1, capsize=2)
plt.errorbar(0., mean, yerr = err_95, marker='x', color='r', elinewidth=0.5, capsize=1)

err_68, err_95 = error_bars(sample)
plt.errorbar(0., mean, yerr = err_68, marker='x', color='r', elinewidth=1, capsize=2)
plt.errorbar(0., mean, yerr = err_95, marker='x', color='r', elinewidth=0.5, capsize=1)


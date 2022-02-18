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
from expected_events_discrete_clustered_optimised import n_ex

# return probability mass function for a Poisson distribution
def poisson_pmf(k, lam):
    return np.power(lam, k) * np.exp(-lam) / factorial(k)

# Range of PBH masses to consider
#m_pbhs = 10**np.arange(0., 1.1, 1)
m_pbhs = 10**np.array([2.])

# Range of numbers of PBHs per cluster to consider
#n_cls = 10**np.arange(7., 2.9, -1.)
n_cls = 10**np.array([5., 6., 7.])

# total number of realisations
n_realisations = 1000


"""Calculate the mean number of expected events"""
for m_pbh in m_pbhs:
    #fig1 = plt.figure()
    #fig1.suptitle('$M_\mathrm{PBH}$ ' + '= {:.0f}'.format(m_pbh) + '$M_\odot$')
    """
    fig2 = plt.figure()
    fig2.suptitle('$M_\mathrm{PBH}$ ' + '= {:.0f}'.format(m_pbh) + '$M_\odot$')
    """
    plt.subplots_adjust(wspace=0.4, hspace=0.4, top=0.85)
   
    for j, n_cl in enumerate(n_cls):
        
        x_min = []
        n_ex_mean = []
        n_ex_poisson = []
        
        x_min_min = 0.05
        i_min = 0
        
        for i in range(n_realisations):
            
            np.random.seed = int(i*m_pbh*n_cl)
            
            filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + str(i)
            d_L = np.loadtxt(filepath + '_dL.txt', delimiter=',')
            v = np.loadtxt(filepath + '_v.txt', delimiter=',')

            # calculate number of expected events
            n_ex_mean.append(n_ex(d_L, v, m_pbh, n_cl, eff=True, poisson=False))
            n_ex_poisson.append(n_ex(d_L, v, m_pbh, n_cl, eff=True, poisson=True))
            
            x_min.append(min(d_L) / 50e3)
            if min(d_L) / 50e3 < x_min_min:
                x_min_min = min(d_L) / 50e3
                i_min = i
                print('i = {:.0f}'.format(i) + ', x={:.5f}'.format(min(d_L)/50e3))

        
        print('N_cl=1e', np.log10(n_cl))
        print('M_PBH=1e', np.log10(m_pbh))

        print('i_min = ', i_min)
        """
        if n_cl > 10**6.5:
            fig3 = plt.figure()
            plt.plot(n_ex_poisson, x_min, 'x')
            plt.xlabel('$N_\mathrm{ex}$')
            plt.ylabel('min$(x)$')
            plt.title('$N_\mathrm{cl} ' + '= 10^{:.0f}$'.format(np.log10(n_cl)) + ', $M_\mathrm{PBH} ' + '= {:.0f}$'.format(m_pbh) + '$M_\odot$')
            plt.yscale('log')
            plt.xscale('log')
            plt.savefig(f'{os.getcwd()}' + '/figures/small_x_test/n_ex_m_pbh=1e{:.0f}_n_cl=1e{:.0f}'.format(np.log10(m_pbh), np.log10(n_cl)) + '.pdf')
            plt.close()
        """
        """
        if n_cl < 10**5.5:
            fig4 = plt.figure()
            axis = np.arange(0, 100, 1)
            plt.hist(n_ex_poisson, density=True, bins=100)
            plt.plot(poisson_pmf(axis, lam=np.mean(n_ex_poisson)), label='Poisson distribution')
            plt.xlabel('$N_\mathrm{ex}$')
            plt.ylabel('$P(N_\mathrm{ex})$')
            plt.legend()
            plt.title('$N_\mathrm{cl} ' + '= 10^{:.0f}$'.format(np.log10(n_cl)) + ', $M_\mathrm{PBH} ' + '= {:.0f}$'.format(m_pbh) + '$M_\odot$')
            plt.savefig(f'{os.getcwd()}' + '/figures/small_x_test/poisson_comparison_n_ex_m_pbh=1e{:.0f}_n_cl=1e{:.0f}'.format(np.log10(m_pbh), np.log10(n_cl)) + '.pdf')
            plt.close()
        """
        """
        ax1 = fig1.add_subplot(2, 2, 1+j)
        ax1.plot(n_ex_mean, x_min, 'x')
        ax1.set_xlabel('$\sum^{N_\mathrm{clusters}}_{i=1} \overline{N}_\mathrm{ex, i}$')
        ax1.set_ylabel('min$(x)$')
        ax1.set_title('$N_\mathrm{cl} ' + '= 10^{:.0f}$'.format(np.log10(n_cl)) + ', $M_\mathrm{PBH} ' + '= {:.0f}$'.format(m_pbh) + '$M_\odot$')
        """
        """
        ax1 = fig1.add_subplot(2, 2, 1+j)
        ax1.plot(n_ex_poisson, x_min, 'x', alpha=0.5)
        ax1.plot(np.mean(n_ex_poisson), np.mean(x_min), 'rx')
        ax1.set_xlabel('$N_\mathrm{ex}$')
        ax1.set_ylabel('min$(x)$')
        ax1.set_title('$N_\mathrm{cl} ' + '= 10^{:.0f}$'.format(np.log10(n_cl)) + ', $M_\mathrm{PBH} ' + '= {:.0f}$'.format(m_pbh) + '$M_\odot$')
        if n_cl < 10**6:
            ax1.set_xlim(0, 50)
            ax1.set_ylim(0, 0.05)
        """
    """
    plt.savefig(f'{os.getcwd()}' + '/figures/small_x_test/n_ex_mean_m_pbh=1e{:.0f}'.format(np.log10(m_pbh)) + '.pdf')
    plt.close()
    """
    plt.savefig(f'{os.getcwd()}' + '/figures/small_x_test/n_ex_m_pbh=1e{:.0f}'.format(np.log10(m_pbh)) + '.pdf')
    plt.close()

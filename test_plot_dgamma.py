#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 15:17:12 2022

@author: ppxmg2
"""

import numpy as np
import matplotlib.pyplot as plt
import halomodel_optimised as hm
import os

# Range of PBH masses to consider
#m_pbhs = 10**np.arange(-1, 3.1, 1)
#m_pbhs = 10**np.arange(0., 0.1, 1)
m_pbhs = np.array([1, 10, 100])

# Uppermost values of the event duration to plot, in days
upper = np.array([300, 1000, 2500])

# Range of numbers of PBHs per cluster to consider
#n_cls = 10**np.arange(3, 7.1, 1.)
n_cls = 10**np.arange(5, 5.1, 1.)

# Number of realisations for each PBH mass and cluster size
n_realisations = 11

# realisation number to start from
min_realisations = 0

# choose a range of realisations to plot
realisations_ncl_1e5_mpbh_1 = np.array([3, 4, 467])
realisations_ncl_1e5_mpbh_10 = np.array([3, 4, 68])
realisations_ncl_1e5_mpbh_100 = np.array([3, 4, 778])


for k, m_pbh in enumerate(m_pbhs):
    
    # import relevant comparison curve for the event duration distribution
    filepath = f'{os.getcwd()}' + '/data_files/event_duration_dist_step=0.50_mpbh={:.2f}'.format(np.log10(m_pbh)) + '.txt'
    t_smooth, dgamma_smooth = np.loadtxt(filepath, delimiter = ',')
        
    for n_cl in n_cls:
        
        fig = plt.figure(figsize=(6, 6))
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.gca().spines['left'].set_visible(False)
        plt.gca().xaxis.set_visible(False)
        plt.gca().yaxis.set_visible(False)
        
        plt.subplots_adjust(wspace=0.5, hspace=0.2, top=0.9)

        setup = hm.Halomodel(m_pbh=m_pbh, n_cl=n_cl)
        
        if np.log10(n_cl) == 5 and np.log10(m_pbh) == 0:
            realisations = realisations_ncl_1e5_mpbh_1
            scale = 1e5
            
        if np.log10(n_cl) == 5 and np.log10(m_pbh) == 1:
            realisations = realisations_ncl_1e5_mpbh_10
            scale = 1e6
            
        if np.log10(n_cl) == 5 and np.log10(m_pbh) == 2:
            realisations = realisations_ncl_1e5_mpbh_100
            scale = 1e7
        
        #for i in range(min_realisations, n_realisations):
        for i, r in enumerate(realisations):
            #ax = fig.add_subplot(n_realisations-min_realisations, 1, i+1-min_realisations)
            ax = fig.add_subplot(3, 1, i+1)
            filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + str(r)
            d_L = np.loadtxt(filepath + '_dL.txt', delimiter=',')
            v = np.loadtxt(filepath + '_v.txt', delimiter=',')
            
            t_hat = 2 * setup.einstein_radius(d_L) / v
            gamma_c = setup.event_rate(d_L, v)
            
            bin_spacing = upper[k] / np.sqrt(len(t_hat))
            
            #scale = 1
                        
            ax.plot(t_smooth, scale * np.array(dgamma_smooth), color='k', linewidth=0.7)
            ax.hist(np.array(t_hat) * 365.25, weights = scale * np.array(gamma_c) / (bin_spacing / 365.25), histtype='step', bins = np.arange(0, upper[k]+bin_spacing/2, bin_spacing))
            ax.set_xlabel('$\hat{t}$ (days)')
            ax.set_ylabel(r'$10^{:.0f}'.format(np.log10(scale)) + '\mathrm{d}\Gamma / \mathrm{d} \hat{t}$ (years)$^{-2}$')
            ax.set_xlim(0, upper[k])
            ax.set_ylim(0, scale * max(dgamma_smooth) * 2.)
            if r == 467:
                ax.set_ylim(0, scale * max(dgamma_smooth) * 3.5)
            if r == 68:
                ax.set_ylim(0, scale * max(dgamma_smooth) * 9.5)
            if r == 778:
                ax.set_ylim(0, scale * max(dgamma_smooth) * 11.)
            #ax.set_title('Realisation ' + str(r))
            
        plt.tight_layout()
        plt.savefig(f'{os.getcwd()}' + '/figures/event_duration_distributions/n_cl=1e{:.0f}_mpbh=1e{:.0f}'.format(np.log10(n_cl), np.log10(m_pbh)) + '.pdf')
        #plt.title('$N_\mathrm{cl}$' + '$ = 10^{:.0f}$'.format(np.log10(n_cl)) + ', $M_\mathrm{PBH} $' + '$= {:.0f} M_\odot$'.format(m_pbh))          
        #plt.savefig(f'{os.getcwd()}' + '/figures/event_duration_distributions/n_cl=1e{:.0f}_mpbh=1e{:.0f}_title'.format(np.log10(n_cl), np.log10(m_pbh)) + '.pdf')


""" Try plotting different N_cl on the same M_PBH plot"""
i = 1
upper = np.array([300, 1000, 2500])
for k, m_pbh in enumerate(m_pbhs):
    fig = plt.figure()
    
    # import relevant comparison curve for the event duration distribution
    filepath = f'{os.getcwd()}' + '/data_files/event_duration_dist_step=0.50_mpbh={:.2f}'.format(np.log10(m_pbh)) + '.txt'
    t_smooth, dgamma_smooth = np.loadtxt(filepath, delimiter = ',')
    plt.plot(t_smooth, dgamma_smooth, color='k', linewidth=0.7)
    
    for n_cl in n_cls:
        
        plt.title('$M_\mathrm{PBH}$' + '$ = {:.0f} M_\odot$'.format(m_pbh))
        setup = hm.Halomodel(m_pbh=m_pbh, n_cl=n_cl)
        
        filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + str(i)
        d_L = np.loadtxt(filepath + '_dL.txt', delimiter=',')
        v = np.loadtxt(filepath + '_v.txt', delimiter=',')
        
        t_hat = 2 * setup.einstein_radius(d_L) / v
        gamma_c = setup.event_rate(d_L, v)
        
        bin_spacing = upper[k] / np.sqrt(len(t_hat))
                    
        plt.hist(np.array(t_hat) * 365.25, weights = np.array(gamma_c) / (bin_spacing / 365.25), histtype='step', bins = np.arange(0, upper[k]+bin_spacing/2, bin_spacing), label='$N_\mathrm{cl} $' + '$= 10^{:.0f}$'.format(np.log10(n_cl)))
        plt.legend()

    plt.xlabel('$\hat{t}$ (days)')
    plt.ylabel(r'$\mathrm{d}\Gamma / \mathrm{d} \hat{t}$ (years)$^{-2}$')
    plt.xlim(0, upper[k])
    plt.ylim(0, max(dgamma_smooth) * 2)
    plt.savefig(f'{os.getcwd()}' + '/figures/event_duration_distributions/n_cl=1e{:.0f}'.format(np.log10(n_cl)) + '.pdf')


""" Check whether the big spikes at short-duration events correspond to small-D_L events"""
n_cl, m_pbh, i, upper = 10000, 10, 3, 300
filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + str(i)
d_L = np.loadtxt(filepath + '_dL.txt', delimiter=',')
v = np.loadtxt(filepath + '_v.txt', delimiter=',')

setup = hm.Halomodel(m_pbh=m_pbh, n_cl=n_cl)
t_hat = 2 * setup.einstein_radius(d_L) / v
gamma_c = setup.event_rate(d_L, v)

plt.figure()
plt.plot(365.25*np.array(t_hat), gamma_c, 'o', alpha=0.3)
plt.yscale('log')

print(np.mean(gamma_c))
print(np.median(gamma_c))

plt.figure()
plt.plot(d_L, gamma_c, 'o', alpha=0.5)
plt.yscale('log')
plt.xscale('log')

print(np.mean(t_hat) * 365.25)

plt.figure()
bin_spacing = 0.01 *  upper / np.sqrt(len(t_hat))
plt.hist(np.array(t_hat) * 365.25, weights = np.array(gamma_c) / (bin_spacing / 365.25), histtype='step', bins = np.arange(0, upper+bin_spacing/2, bin_spacing), label='$N_\mathrm{cl} $' + '$= 10^{:.0f}$'.format(np.log10(n_cl)))

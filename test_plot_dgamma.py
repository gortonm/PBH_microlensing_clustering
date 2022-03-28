# Code to plot Fig. 1 (event duration distributions)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import halomodel_optimised as hm
import os

# Specify the plot style
mpl.rcParams.update({'font.size': 14,'font.family':'serif'})
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

# Range of PBH masses to consider
m_pbhs = np.array([1, 10])

# Uppermost values of the event duration to plot, in days
upper = np.array([300, 1000, 2500])

# Range of numbers of PBHs per cluster to consider
#n_cls = 10**np.arange(3, 7.1, 1.)
n_cls = 10**np.arange(6, 6.1, 1.)    # for Fig. 1, only show N_cl = 1e5

# choose realisations to plot (these depend on the PBH mass)
realisations_ncl_1e5_mpbh_1 = np.array([0, 1, 607])
#realisations_ncl_1e5_mpbh_1 = np.array([116, 125, 135])
realisations_ncl_1e5_mpbh_10 = np.array([0, 1, 847])
realisations = np.array([0,1,2])
scale=1

for k, m_pbh in enumerate(m_pbhs):
    
    # import relevant comparison curve for the event duration distribution
    filepath = f'{os.getcwd()}' + '/data_files/event_duration_dist_step=0.50_mpbh={:.2f}'.format(np.log10(m_pbh)) + '.txt'
    t_smooth, dgamma_smooth = np.loadtxt(filepath, delimiter = ',')
        
    for n_cl in n_cls:
        
        fig = plt.figure(figsize=(4.5, 7))
        gs = fig.add_gridspec(3, ncols=1, hspace=0)
        axes = gs.subplots(sharex=True)
        
        setup = hm.Halomodel(m_pbh=m_pbh, n_cl=n_cl)
        
        # scale = factor to divide y-axis values by to show more clearly on plot
        if np.log10(n_cl) == 6 and np.log10(m_pbh) == 0:
            realisations = realisations_ncl_1e5_mpbh_1
            scale = 1e5
            
        if np.log10(n_cl) == 6 and np.log10(m_pbh) == 1:
            realisations = realisations_ncl_1e5_mpbh_10
            scale = 1e6
                
        for i, r in enumerate(realisations):
            
            ax = axes[i]
            
            filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + str(r)
            d_L = np.loadtxt(filepath + 'updateRE_dL.txt', delimiter=',')
            v = np.loadtxt(filepath + 'updateRE_v.txt', delimiter=',')
            
            t_hat = 2 * setup.einstein_radius(d_L) / v
            gamma_c = setup.event_rate(d_L, v)
            
            bin_spacing = upper[k] / np.sqrt(len(t_hat))
                                    
            ax.set_xticks(np.arange(0, upper[k], 6))
            ax.plot(t_smooth, scale * np.array(dgamma_smooth), color='k', linewidth=0.7)
            ax.hist(np.array(t_hat) * 365.25, weights = scale * np.array(gamma_c) / (bin_spacing / 365.25), histtype='step', bins = np.arange(0, upper[k]+bin_spacing/2, bin_spacing))
            ax.set_xlim(0, upper[k])
            ax.set_ylim(0, scale * max(dgamma_smooth) * 1.5)
            ax.set_xticks(np.linspace(0, upper[k], 6))
            
            # Values of r correspond to specific realisations with a large spike at small event durations
            # For these cases, increase the maximum value of the y-axis
            if r == 607:
                ax.set_ylim(0, scale * max(dgamma_smooth) * 9.5)
            if r == 847:
                ax.set_ylim(0, scale * max(dgamma_smooth) * 16)

        # add an invisible axis, for the axis labels that apply for the whole figure
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
            
        plt.xlabel('$\hat{t}$ (days)', fontsize='16')
        plt.ylabel(r'$10^{:.0f}'.format(np.log10(scale)) + '\mathrm{d}\Gamma / \mathrm{d} \hat{t}$ (years)$^{-2}$', fontsize='16', labelpad=10)

        plt.tight_layout()
        plt.savefig(f'{os.getcwd()}' + '/figures/event_duration_distributions/n_cl=1e{:.0f}_mpbh=1e{:.0f}'.format(np.log10(n_cl), np.log10(m_pbh)) + '_smallestx.pdf')
# Code to plot Fig. 1 of paper (differential event rate)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

# Specify the plot style
mpl.rcParams.update({'font.size': 14, 'font.family': 'serif'})
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


# If True, use use EROS-2 efficiency function and exposure
use_eff_EROS2 = True

# If True, use 'toy' survey efficiency function and exposure
use_eff_toy = True

# if True, set PBH cluster radius to 10 pc
set_rcl_10 = True

append = ""
if set_rcl_10:
    append += "_rcl10"
if use_eff_EROS2:
    append += "_EROS2"

# Fraction of dark matter in PBHs
f_pbh = 1

# PBH masses (in solar masses)
m_pbhs = np.array([1, 10])

# Uppermost values of the event duration to plot, in days
upper = np.array([300, 1000])

# Range of numbers of PBHs per cluster to consider
n_cl = 1e6    # for Fig. 1, only show N_cl = 1e6

# Realisations with a cluster close to the observer
r_small_x_m_pbh_1, r_small_x_m_pbh_10 = 93, 2239
# Realisations to plot (these depend on the PBH mass)
realisations_ncl_1e6_mpbh_1 = np.array([3, 4, r_small_x_m_pbh_1])
realisations_ncl_1e6_mpbh_10 = np.array([2, 3, r_small_x_m_pbh_10])

for k, m_pbh in enumerate(m_pbhs):

    # import relevant comparison curve for the event duration distribution
    filepath = f'{os.getcwd()}' + \
        '/data_files/event_duration_dist_step=0.50_mpbh={:.2f}'.format(
            np.log10(m_pbh)) + '.txt'
    t_smooth, dgamma_smooth = np.loadtxt(filepath, delimiter=',')

    fig = plt.figure(figsize=(4.5, 7))
    gs = fig.add_gridspec(3, ncols=1, hspace=0)
    axes = gs.subplots(sharex=True)

    # scale = factor to divide y-axis values by to show more clearly on plot
    if m_pbh == 1:
        realisations = realisations_ncl_1e6_mpbh_1
        scale = 1e5

    if m_pbh == 10:
        realisations = realisations_ncl_1e6_mpbh_10
        scale = 1e6

    for i, r in enumerate(realisations):

        # Load data from csv file
        filename = f'{os.getcwd()}' + '/simulated_data/log10_Ncl/{0:.2f}'.format(np.log10(n_cl)) + \
            '/log10_mpbh/{0:.2f}/'.format(np.log10(m_pbh)) + \
            "cluster_data_" + str(r) + append + ".csv"

        data = np.genfromtxt(filename, delimiter=',')
        data_transpose = np.transpose(data)

        dL_values = data_transpose[0]
        v_values = data_transpose[1]
        # distance from cluster centre from line of sight
        d_values = data_transpose[2]
        event_rate_values = data_transpose[3]
        t_hat_values_corrected = data_transpose[4]
        f_cl = data_transpose[5]

        # For 'corrected' subsample of event rates, include the area correction factor
        event_rate_values_corrected = f_cl * event_rate_values

        ax = axes[i]
        bin_spacing = upper[k] / np.sqrt(len(t_hat_values_corrected))

        ax.set_xticks(np.arange(0, upper[k], 6))
        ax.plot(t_smooth, scale * np.array(dgamma_smooth),
                color='k', linewidth=0.7)
        ax.hist(np.array(t_hat_values_corrected) * 365.25, weights=scale * np.array(event_rate_values_corrected) /
                (bin_spacing / 365.25), histtype='step', bins=np.arange(0, upper[k]+bin_spacing/2, bin_spacing))
        ax.set_xlim(0, upper[k])
        ax.set_ylim(0, scale * max(dgamma_smooth) * 1.5)
        ax.set_xticks(np.linspace(0, upper[k], 6))

        # Values of r correspond to specific realisations with a large spike at small event durations
        # For these cases, increase the maximum value of the y-axis
        if r == r_small_x_m_pbh_1:
            ax.set_ylim(0, scale * max(dgamma_smooth) * 8.5)
        if r == r_small_x_m_pbh_10:
            ax.set_ylim(0, scale * max(dgamma_smooth) * 12)

    # Add an invisible axis, for the axis labels that apply for the whole figure
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False,
                    bottom=False, left=False, right=False)
    plt.xlabel('$\hat{t}$ (days)', fontsize='16')
    plt.ylabel(r'$10^{:.0f}'.format(np.log10(
        scale)) + '\mathrm{d}\Gamma / \mathrm{d} \hat{t}$ (years)$^{-2}$', fontsize='16', labelpad=10)
    plt.tight_layout()
    plt.savefig(f'{os.getcwd()}' + '/figures/differential_event_rate_Ncl=1e{:.0f}_Mpbh={:.0f}'.format(
        np.log10(n_cl), m_pbh) + 'Msun' + '.pdf')

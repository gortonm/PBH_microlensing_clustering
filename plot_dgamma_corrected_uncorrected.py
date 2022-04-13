# Code to plot Fig. 1 (event duration distributions)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from generate_results_updated_pcl_v2 import r_cone

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


# Boolean controlling whether to use EROS-2 efficiency function and exposure
EROS2_eff = True

# Boolean controlling whether to set PBH cluster radius to 10 pc
set_rcl_10 = True

append = ""
if set_rcl_10:
    append += "_rcl10_"
if EROS2_eff:
    append += "_EROS2_"

f_pbh = 1
x_max = 1
d_s = 50e3

# Range of PBH masses to consider
m_pbhs = np.array([1])

# Uppermost values of the event duration to plot, in days
upper = np.array([300, 1000])

# Range of numbers of PBHs per cluster to consider
n_cl = 1e6    # for Fig. 1, only show N_cl = 1e5

r_small_m_pbh_1, r_small_m_pbh_10 = 93, 2239
# choose realisations to plot (these depend on the PBH mass)
#realisations_ncl_1e6_mpbh_1 = np.arange(0, 5, 1)
realisations_ncl_1e6_mpbh_1 = np.array([93, 684, 1034, 3193, 3333, 5029, 7882, 8354, 8788])
#realisations_ncl_1e6_mpbh_1 = np.array([644, 724, 843, 1061, 1949, 2496, 2560, 2977])
realisations_ncl_1e6_mpbh_1 = np.array([4552, 4795, 4991, 5032, 5380, 6003, 8140, 9733])
#realisations_ncl_1e6_mpbh_1 = np.append(realisations_ncl_1e6_mpbh_1, r_small_m_pbh_1)

realisations_ncl_1e6_mpbh_10 = np.arange(0, 5, 1)
realisations_ncl_1e6_mpbh_10 = np.arange(93, 684, 1034)

#np.append(realisations_ncl_1e6_mpbh_10, r_small_m_pbh_10)
scale=1


realisations_ncl_1e6_mpbh_1 = np.array([1, 2, 3, 4])
realisations_ncl_1e6_mpbh_1 = np.array([724, 843, 1949])
realisations_ncl_1e6_mpbh_1 = np.array([93, 684, 1034])


for k, m_pbh in enumerate(m_pbhs):
    
    # import relevant comparison curve for the event duration distribution
    filepath = f'{os.getcwd()}' + '/data_files/event_duration_dist_step=0.50_mpbh={:.2f}'.format(np.log10(m_pbh)) + '.txt'
    t_smooth, dgamma_smooth = np.loadtxt(filepath, delimiter = ',')
        
    fig = plt.figure(figsize=(4.5, 7))
    gs = fig.add_gridspec(len(realisations_ncl_1e6_mpbh_1), ncols=1, hspace=0)
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
        filename = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + "cluster_data_update_pcl_" + str(r) + append + ".csv"    
        
        data = np.genfromtxt(filename, delimiter=',')
        data_transpose = np.transpose(data)
        
        dL_values_corrected_import = data_transpose[0]
        v_values_corrected_import = data_transpose[1]
        d_values_import = data_transpose[2]
        event_rate_values_import = data_transpose[3]
        t_hat_values_corrected_import = data_transpose[4]
        f_cl_import = data_transpose[5]
        
        # Restrict to clusters with a fractional line of sight distance smaller than x_max
        # Set x_max = 1 to include all clusters
        dL_values_corrected = dL_values_corrected_import[dL_values_corrected_import < x_max * d_s]
        v_values_corrected = v_values_corrected_import[dL_values_corrected_import < x_max * d_s]
        f_cl = f_cl_import[dL_values_corrected_import < x_max * d_s]
        d_values = d_values_import[dL_values_corrected_import < x_max * d_s]
        event_rate_values = event_rate_values_import[dL_values_corrected_import < x_max * d_s]
        t_hat_values_corrected = t_hat_values_corrected_import[dL_values_corrected_import < x_max * d_s]
    
        # comment out to produce Figure 1
        
        # 'Uncorrected' subsample: include only clusters with cluster centres located within the microlensing cone
        dL_values_uncorrected = dL_values_corrected[ d_values < r_cone(dL_values_corrected) ]
        v_values_uncorrected = v_values_corrected[ d_values < r_cone(dL_values_corrected) ]
        event_rate_values_uncorrected = event_rate_values[ d_values < r_cone(dL_values_corrected) ]
        t_hat_values_uncorrected = t_hat_values_corrected[ d_values < r_cone(dL_values_corrected) ]
        
        # For 'corrected' subsample of event rates, include the area correction factor
        event_rate_values_corrected = f_cl * event_rate_values
        
        ax = axes[i]
        
        bin_spacing = upper[k] / np.sqrt(len(t_hat_values_corrected))
                                
        ax.set_xticks(np.arange(0, upper[k], 6))
        ax.plot(t_smooth, scale * np.array(dgamma_smooth), color='k', linewidth=0.7)
        ax.hist(np.array(t_hat_values_corrected) * 365.25, weights = scale * np.array(event_rate_values_corrected) / (bin_spacing / 365.25), color='tab:blue', histtype='step', bins = np.arange(0, upper[k]+bin_spacing/2, bin_spacing))
        ax.hist(np.array(t_hat_values_uncorrected) * 365.25, weights = scale * np.array(event_rate_values_uncorrected) / (bin_spacing / 365.25), color='tab:orange', histtype='step', linewidth = 0.5, bins = np.arange(0, upper[k]+bin_spacing/2, bin_spacing))   
        ax.set_xlim(0, upper[k])
        ax.set_ylim(0, scale * max(dgamma_smooth) * 18)
        #ax.set_ylim(0, scale * max(dgamma_smooth) * 2)
        ax.set_xticks(np.linspace(0, upper[k], 6))
        
        if r== 724:
            ax.set_ylim(0, scale * max(dgamma_smooth) * 5)
            
        if r == 684:
            ax.set_ylim(0, scale * max(dgamma_smooth) * 14)

        
        # Values of r correspond to specific realisations with a large spike at small event durations
        # For these cases, increase the maximum value of the y-axis
        if r == r_small_m_pbh_1:
            ax.set_ylim(0, scale * max(dgamma_smooth) * 14)
        if r == r_small_m_pbh_10:
            ax.set_ylim(0, scale * max(dgamma_smooth) * 12)

    # add an invisible axis, for the axis labels that apply for the whole figure
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
            
    plt.xlabel('$\hat{t}$ (days)', fontsize='16')
    plt.ylabel(r'$10^{:.0f}'.format(np.log10(scale)) + '\mathrm{d}\Gamma / \mathrm{d} \hat{t}$ (years)$^{-2}$', fontsize='16', labelpad=10)

    # Add legend outside of figure with descriptions of the curves
    line_corrected = mpl.lines.Line2D([], [], color='tab:blue', label='Cluster centres in extended region \n (with area correction)')
    line_unccorrected = mpl.lines.Line2D([], [], linewidth=0.5, color='tab:orange', label='Cluster centres in microlensing cone \n (no area correction)')
    #plt.legend(handles=[line_corrected, line_unccorrected], fontsize='small', bbox_to_anchor=(1,1), loc="upper left")

    plt.tight_layout()
    #plt.savefig(f'{os.getcwd()}' + '/figures/n_cl=1e{:.0f}_mpbh=1e{:.0f}'.format(np.log10(n_cl), np.log10(m_pbh)) + '_smallx_only_in_corrected.pdf')
    plt.savefig(f'{os.getcwd()}' + '/figures/n_cl=1e{:.0f}_mpbh=1e{:.0f}'.format(np.log10(n_cl), np.log10(m_pbh)) + '_smallx_in_uncorrected.pdf')
    #plt.savefig(f'{os.getcwd()}' + '/figures/n_cl=1e{:.0f}_mpbh=1e{:.0f}'.format(np.log10(n_cl), np.log10(m_pbh)) + '.pdf')
    
    
    
"""Plot a single realisation"""
m_pbh = 1
r = 1949


# Load data from csv file
filename = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + "cluster_data_update_pcl_" + str(r) + append + ".csv"    

data = np.genfromtxt(filename, delimiter=',')
data_transpose = np.transpose(data)

dL_values_corrected_import = data_transpose[0]
v_values_corrected_import = data_transpose[1]
d_values_import = data_transpose[2]
event_rate_values_import = data_transpose[3]
t_hat_values_corrected_import = data_transpose[4]
f_cl_import = data_transpose[5]

# Restrict to clusters with a fractional line of sight distance smaller than x_max
# Set x_max = 1 to include all clusters
dL_values_corrected = dL_values_corrected_import[dL_values_corrected_import < x_max * d_s]
v_values_corrected = v_values_corrected_import[dL_values_corrected_import < x_max * d_s]
f_cl = f_cl_import[dL_values_corrected_import < x_max * d_s]
d_values = d_values_import[dL_values_corrected_import < x_max * d_s]
event_rate_values = event_rate_values_import[dL_values_corrected_import < x_max * d_s]
t_hat_values_corrected = t_hat_values_corrected_import[dL_values_corrected_import < x_max * d_s]

# comment out to produce Figure 1

# 'Uncorrected' subsample: include only clusters with cluster centres located within the microlensing cone
dL_values_uncorrected = dL_values_corrected[ d_values < r_cone(dL_values_corrected) ]
v_values_uncorrected = v_values_corrected[ d_values < r_cone(dL_values_corrected) ]
event_rate_values_uncorrected = event_rate_values[ d_values < r_cone(dL_values_corrected) ]
t_hat_values_uncorrected = t_hat_values_corrected[ d_values < r_cone(dL_values_corrected) ]

# For 'corrected' subsample of event rates, include the area correction factor
event_rate_values_corrected = f_cl * event_rate_values

plt.figure()
ax = plt.gca()

bin_spacing = 300 / (365.25 * np.sqrt(len(t_hat_values_corrected)))
                        
ax.plot(np.array(t_smooth)/365.25, scale * np.array(dgamma_smooth), color='k', linewidth=0.7)
ax.hist(np.array(t_hat_values_corrected), weights = scale * np.array(event_rate_values_corrected) / (bin_spacing), color='tab:blue', histtype='step', bins = np.arange(0, (300/365.25)+bin_spacing/2, bin_spacing))
ax.hist(np.array(t_hat_values_uncorrected), weights = scale * np.array(event_rate_values_uncorrected) / (bin_spacing), color='tab:orange', histtype='step', linewidth = 0.5, bins = np.arange(0, (300/365.25)+bin_spacing/2, bin_spacing))        
ax.set_ylim(0, scale * max(dgamma_smooth) * 2)
plt.xlim(0, 300/365.25)
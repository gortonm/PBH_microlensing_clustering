# Code to plot Fig. 2 (P(N_obs) against N_obs)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from scipy.special import factorial


# Specify the plot style
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


def poisson_pmf(k, lam):
    """
    Probability mass function for a Poisson distribution

    Parameters
    ----------
    k : Integer
        Number of events.
    lam : Float
        Mean number of events.

    Returns
    -------
    Float
        Probability of k events, for a Poisson-distributed variable with a mean
        lam events.

    """
    return np.power(lam, k) * np.exp(-lam) / factorial(k)


new_RE = False
seed_i = True
no_save_perfect = True
extended_range_mpbh = True

if seed_i:
    a = '_seed_i'
else:
    a = ''
    
if extended_range_mpbh:
    a = 'extended_range_MPBH'
    
if no_save_perfect:
    a += '_nosaveperfect'
    
if new_RE:
    append = 'newRE' + a
    
else:
    append = 'oldRE' + a

# PBH masses
#m_pbhs = np.array([1, 100])
m_pbhs = np.array([1])

# Nimbers of PBHs per cluster to consider
#n_cls = 10**np.arange(4, 8.1, 1.)
n_cls = 10**np.arange(7, 8.1, 1.)   # uncomment if only using N_cl = 1e6, 1e7, 1e8

# Number of realisations for each PBH mass and number of PBHs per cluster
n_realisations = 10000

# Colours to display for different N_cl values
#colours = ['b', 'm', 'darkorange', 'r', 'saddlebrown']
colours = ['darkorange', 'r', 'saddlebrown']    # uncomment if only using N_cl = 1e6, 1e7, 1e8

# Uppermost value of N_obs to display in plot
n_ex_upper = [150, 12]


for j, m_pbh in enumerate(m_pbhs):
    
    # Mean value of the number of expected events, from smoothly-distributed PBHs of mass m_pbh
    # This is calculated by numerically integrating the mean number of events for the smooth case
    if np.log10(m_pbh) == 0:
        mean_smooth = 24.6
        mean_smooth_perfect = 60.1
        
    if np.log10(m_pbh) == 2:
        mean_smooth = 0.278
        mean_smooth_perfect = 6.0        
    
    # Set up figure and bins for histogram
    bins = np.arange(0, 1000, 1)
    plt.figure(figsize=(5, 4))
    
    # Add Poisson probability mass function for mean value from smooth case
    plt.plot(poisson_pmf(np.arange(0, 101, 1.), lam = round(mean_smooth)), color='k', linestyle='dashed')
    plt.plot(poisson_pmf(np.arange(0, 101, 1.), lam = round(mean_smooth_perfect)), color='k', linestyle='dotted')
    
    print(round(mean_smooth))
    round(round(mean_smooth_perfect))
    
    # maximum values of i such that m_pbh * n_cl * i < 2**32 - 1
    n = (10000, 429, 42)
    
    for i, n_cl in enumerate(n_cls):        
        filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh))
        # Load file for number of events in all realisations, for a given PBH mass and number of PBHs per cluster
        n_ex_EROS_efficiency = np.loadtxt(filepath + 'n_ex_EROS_2_fpbh=1.000_1e4samples' + append + '.txt')
        n_ex_perfect_efficiency = np.loadtxt(filepath + 'n_ex_perfect_fpbh=1.000_1e4samples' + append + '.txt')
        
        # replace the first values of the outputs with the result obtained if using a seed n_cl * m_pbh * i
        if seed_i:
            n_ex_EROS_efficiency = np.loadtxt(filepath + 'n_ex_EROS_2_fpbh=1.000_1e4samples' + 'oldRE' + '.txt')
            n_ex_perfect_efficiency = np.loadtxt(filepath + 'n_ex_perfect_fpbh=1.000_1e4samples' + 'oldRE' + '.txt')
            
            n_ex_EROS_efficiency_seed_i = np.loadtxt(filepath + 'n_ex_EROS_2_fpbh=1.000_1e4samples' + append + '.txt')
            #n_ex_perfect_efficiency_seed_i = np.loadtxt(filepath + 'n_ex_perfect_fpbh=1.000_1e4samples' + append + '.txt')
            
            k = 0
            while k < n[i]:
                n_ex_EROS_efficiency[k] = n_ex_EROS_efficiency_seed_i[k]
                k += 1 
        """
        else:
            n_ex_EROS_efficiency = np.loadtxt(filepath + 'n_ex_EROS_2_fpbh=1.000_1e4samples' + append + '.txt')
            n_ex_perfect_efficiency = np.loadtxt(filepath + 'n_ex_perfect_fpbh=1.000_1e4samples' + append + '.txt')
        """
        plt.hist(n_ex_EROS_efficiency, bins=bins, density=True, color=colours[i], histtype='step', label='$N_\mathrm{cl} = $'+'$10^{:.0f}$'.format(np.log10(n_cl)))
        plt.hist(n_ex_perfect_efficiency, bins=bins, density=True, color=colours[i], histtype='step', linestyle='dotted')
                 
    plt.xlim(0, n_ex_upper[j])
    plt.xlabel('$N_\mathrm{obs}$')
    plt.ylabel('$P[N_\mathrm{obs}(f=1)]$')
    plt.tight_layout()
    plt.legend()
            
    plt.savefig(f'{os.getcwd()}' + '/figures/P(N_obs)_Mpbh={:.0f}'.format((m_pbh)) + 'Msun_1e4samples' + append + '.pdf')
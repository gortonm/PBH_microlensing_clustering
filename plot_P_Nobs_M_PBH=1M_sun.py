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


# Boolean controlling whether to use EROS-2 efficiency function and exposure
EROS2_eff = False

# Boolean controlling whether to set PBH cluster radius to 10 pc
set_rcl_10 = True

append = ""
if set_rcl_10:
    append += "_rcl10_"
if EROS2_eff:
    append += "_EROS2_"

f_pbh = 1

# PBH mass
m_pbh = 1

if m_pbh == 1:
    mean_smooth = 24.6
    mean_smooth_perfect = 60.1

# Nimbers of PBHs per cluster to consider
n_cls = 10**np.arange(6, 8.1, 1.)

# Number of realisations for each PBH mass and number of PBHs per cluster
n_realisations = 10000

# Colours to display for different N_cl values
colours = ['darkorange', 'r', 'saddlebrown']

# Uppermost value of N_obs to display in plot
n_ex_upper = 100
        
# Set up figure and bins for histogram
bins = np.arange(0, 1000, 1)
plt.figure(figsize=(5.5, 4.5))


# Add Poisson probability mass function for mean value from smooth case    
for i, n_cl in enumerate(n_cls):     
    filepath = (
        f"{os.getcwd()}"
        + "/simulated_data_constraints/N_cl/{0:.2f}".format(np.log10(n_cl))
        + "/M_PBH/{0:.2f}/".format(np.log10(m_pbh))
    )
    
    # Load file for number of events in all realisations, for a given PBH mass and number of PBHs per cluster
    filename_Nobs_corrected = (
        filepath
        + "n_ex_corrected_fpbh={0:.4f}".format(f_pbh)
        + append
        + "_nsamp={0:.1e}".format(n_realisations)
        + ".txt"
    )
    
    filename_Nobs_perfect = (
        filepath
        + "n_ex_perfect_fpbh={0:.4f}".format(f_pbh)
        + append
        + "_nsamp={0:.1e}".format(n_realisations)
        + ".txt"
    )
    
    
    filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh))
    n_ex_EROS_efficiency = np.loadtxt(filepath + filename_Nobs_corrected)
    n_ex_perfect_efficiency = np.loadtxt(filepath + filename_Nobs_perfect)
    plt.hist(n_ex_EROS_efficiency, bins=bins, density=True, color=colours[i], histtype='step', label='$N_\mathrm{cl} = $'+'$10^{:.0f}$'.format(np.log10(n_cl)))
    plt.hist(n_ex_perfect_efficiency, bins=bins, density=True, color=colours[i], histtype='step', linestyle='dotted')

plt.plot(poisson_pmf(np.arange(0, 101, 1.), lam = mean_smooth), color='k', linestyle='dashed')
plt.plot(poisson_pmf(np.arange(0, 101, 1.), lam = mean_smooth_perfect), color='k', linestyle='dotted')

plt.xlim(0, n_ex_upper)
plt.xlabel('$N_\mathrm{obs}$')
plt.ylabel('$P[N_\mathrm{obs}(f_\mathrm{PBH} = 1)]$')
plt.tight_layout()
plt.legend()
plt.savefig('P(N_obs)_Mpbh={:.0f}'.format((m_pbh)) + 'Msun.pdf')
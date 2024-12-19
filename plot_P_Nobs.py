# Code to plot Fig. 2 of paper (P(N_obs) against N_obs)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from scipy.special import factorial


# Specify the plot style
mpl.rcParams.update({'font.size': 16, 'font.family': 'serif'})
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
    Probability mass function for a Poisson distribution.

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


# If True, use use EROS-2 efficiency function and exposure
use_eff_EROS2 = False

# If True, use 'toy' survey efficiency function and exposure
use_eff_toy = True

# if True, set PBH cluster radius to 10 pc
set_rcl_10 = True

append = ""
if set_rcl_10:
    append += "_rcl10"

m_pbh = 1000   # PBH mass (in solar masses)

# Number of realisations for each PBH mass and number of PBHs per cluster
n_realisations = 10000

# Nimbers of PBHs per cluster to consider
if m_pbh == 1:
    # Numbers of PBHs per cluster to consider
    n_cls = 10**np.arange(6, 8.1, 1.)
    n_ex_upper = 150  # Uppermost value of N_obs to display in plot
    colours = ['darkorange', 'r', 'saddlebrown']

    # Mean value of the number of expected events, from smoothly-distributed
    # PBHs of mass m_pbh, calculated by numerically integrating the mean number
    # of events for the smooth case
    mean_smooth = 24.6
    mean_smooth_perfect = 60.1

elif m_pbh == 1000:
    # Numbers of PBHs per cluster to consider
    n_cls = 10**np.arange(3, 5.1, 1.)
    n_ex_upper = 100  # Uppermost value of N_obs to display in plot
    colours = ['limegreen', 'b', 'm']

    # Mean value of the number of expected events, from smoothly-distributed
    # PBHs of mass m_pbh, calculated by numerically integrating the mean number
    # of events for the smooth case
    mean_smooth = 39.5

# Set up figure and bins for histogram
bins = np.arange(0, 1000, 1)
plt.figure(figsize=(5, 4))

# Add Poisson probability mass function for mean value from smooth case
plt.plot(poisson_pmf(np.arange(0, 101, 1.), lam=round(
    mean_smooth)), color='k', linestyle='dashed')
print("Mean number of events (smooth DM) = %s" % round(mean_smooth))

if m_pbh == 1:
    plt.plot(poisson_pmf(np.arange(0, 101, 1.), lam=round(
        mean_smooth_perfect)), color='k', linestyle='dotted')
    print("Mean number of events (smooth DM, perfect efficiency) = %s"
          % round(mean_smooth_perfect))

for i, n_cl in enumerate(n_cls):
    # Load file for number of events in all realisations, for a given PBH mass
    # and number of PBHs per cluster
    filepath = f'{os.getcwd()}' + '/simulated_data/log10_Ncl/{0:.2f}'.format(
        np.log10(n_cl)) + '/log10_mpbh/{0:.2f}/'.format(np.log10(m_pbh))

    if m_pbh == 1:
        n_ex_EROS_efficiency = np.loadtxt(
            filepath
            + 'n_ex_fpbh=1.0000'
            + append
            + "_EROS2_nsamp={0:.1e}".format(n_realisations)
            + '.txt')
        n_ex_perfect_efficiency = np.loadtxt(
            filepath
            + 'n_ex_perfect_fpbh=1.0000'
            + append
            + "_EROS2_nsamp={0:.1e}".format(n_realisations)
            + '.txt')
        plt.hist(n_ex_EROS_efficiency,
                 bins=bins,
                 density=True,
                 color=colours[i],
                 histtype='step',
                 label=r'$N_\mathrm{cl}=$'+'$10^{:.0f}$'.format(np.log10(n_cl)))
        plt.hist(n_ex_perfect_efficiency,
                 bins=bins,
                 density=True,
                 color=colours[i],
                 histtype='step',
                 linestyle='dotted')

    elif m_pbh == 1000:
        n_ex_toy_efficiency = np.loadtxt(
            filepath
            + 'n_ex_fpbh=1.0000'
            + append
            + "_nsamp={0:.1e}".format(n_realisations)
            + '.txt')
        plt.hist(n_ex_toy_efficiency,
                 bins=bins,
                 density=True,
                 color=colours[i],
                 histtype='step',
                 linestyle="solid",
                 label=r'$N_\mathrm{cl}=$'+'$10^{:.0f}$'.format(np.log10(n_cl)))

plt.xlim(0, n_ex_upper)
plt.xlabel(r'$N_\mathrm{obs}$')
plt.ylabel(r'$P[N_\mathrm{obs}(f=1)]$')
plt.tight_layout()
plt.legend()

plt.savefig(f'{os.getcwd()}'
            + '/figures/P(N_obs)_Mpbh={:.0f}'.format((m_pbh))
            + 'Msun_1e4samples' + '.pdf')

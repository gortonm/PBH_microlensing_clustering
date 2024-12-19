#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 14:03:30 2022

@author: ppxmg2
"""

import numpy as np
import os
import matplotlib.pyplot as plt

# If True, use use EROS-2 efficiency function and exposure
use_eff_EROS2 = False

# If True, use 'toy' survey efficiency function and exposure
use_eff_toy = True

# if True, set PBH cluster radius to 10 pc
set_rcl_10 = True

append = ""
if set_rcl_10:
    append += "_rcl10"
if use_eff_EROS2:
    append += "_EROS2"

# Set number of PBHs per cluster and PBH mass (in solar masses)
n_cl, m_pbh = 1e3, 1e3
# Number of realisations in Monte Carlo simulations
n_realisations = 300000


def prob_n_obs(n_obs_data, n_obs):
    """
    Find probability corresponding to a given value of the number of observed
    events.

    Parameters
    ----------
    n_obs_data : Array-like
        Simulated data for the number of observed microlensing events in each
        realisation of the Monte Carlo simulations.
    n_obs : Float
        Number of microlensing events.

    Returns
    -------
    Float
        Probability of observing a given number of microlensing events.

    """

    # Create histogram to approximate full probability distribution function
    bins = np.arange(-0.5, max(n_obs_data)+0.51, 1)

    probability, bin_edges = np.histogram(n_obs_data, bins=bins, density=True)

    # Since bins are centred at integer values, and the first bin is centred at
    # zero, the index of the bin matches that of the integer value of the
    # number of events
    if n_obs < len(probability):
        return probability[n_obs]
    else:
        return 0


# Fraction of dark matter in PBHs
f_pbhs = np.arange(0.0945, 0.0976, 0.0001)

for f_pbh in f_pbhs:

    n_obs_data = np.loadtxt(
        f"{os.getcwd()}"
        + "/simulated_data/log10_Ncl/{0:.2f}".format(np.log10(n_cl))
        + "/log10_mpbh/{0:.2f}/".format(np.log10(m_pbh))
        + 'n_ex_fpbh={0:.4f}'.format(f_pbh)
        + append
        + '_nsamp={0:.1e}'.format(n_realisations) + '.txt')

    print('f_PBH = {:.4f}, P(N_obs=0) = {:.4f}'.format(
        f_pbh, prob_n_obs(n_obs_data, 0)))

# %% Plot some of the probability distributions of the number of events

# Fraction of dark matter in PBHs
f_pbhs = np.arange(0.095, 0.0976, 0.001)

fig, ax = plt.subplots(figsize=(6, 5))
n_obs_values = np.arange(0, 11, 1)

for f_pbh in f_pbhs:

    n_obs_data = np.loadtxt(
        f"{os.getcwd()}"
        + "/simulated_data/log10_Ncl/{0:.2f}".format(np.log10(n_cl))
        + "/log10_mpbh/{0:.2f}/".format(np.log10(m_pbh))
        + 'n_ex_fpbh={0:.4f}'.format(f_pbh)
        + append
        + '_nsamp={0:.1e}'.format(n_realisations) + '.txt')

    print('f_PBH = {:.4f}, P(N_obs=0) = {:.4f}'.format(
        f_pbh, prob_n_obs(n_obs_data, 0)))
    prob_n_obs_values = [prob_n_obs(n_obs_data, n_obs)
                         for n_obs in n_obs_values]
    ax.plot(n_obs_values, prob_n_obs_values,
            label=r"$f_{\rm PBH} " + "= {:.4f}$".format(f_pbh))

ax.legend()
plt.xlim(0, max(n_obs_values))
plt.xlabel(r'$N_\mathrm{obs}$')
plt.ylabel(r'$P[N_\mathrm{obs}]$')
plt.xticks(n_obs_values)
plt.tight_layout()

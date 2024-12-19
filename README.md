# PBH_microlensing_clustering


Code for calculating microlensing constraints for clustered primordial black hole (PBH) dark matter (DM), associated with the paper "*Effect of clustering on primordial black hole microlensing constraints*"  (see https://arxiv.org/abs/2203.04209). We use Monte Carlo simulations to estimate the probability distribution of the number of expected events, from which a 95% confidence limit on the abundance of PBHs can be estimated.


# Contents
The code, plots etc. are arranged as follows

# Motivation

Stellar microlensing is the temporary magnification of a star that occurs when a compact object passes close to the line of sight to the star. Various microlensing surveys have placed tight constraints on primordial black holes (PBHs), which if taken at face value exclude PBHs of mass $10^{-11} M_\odot \lesssim M_{\rm PBH} \lesssim 10^4 M_\odot$ from making up all of the dark matter (DM). These constraints have been obtained assuming that the DM is smoothly distributed. This code investigates the effect of PBH clustering on stellar microlensing constraints, for PBHs formed from the collapse of large gaussian density fluctuations.

## Folders

* `/data_files` Includes files required to calculate constraints for the case of smoothly-distributed DM (which has been previously considered in the literature and provides a comparison) and the differential event rate for 
smoothly-distributed DM.
* `/figures` Includes figures produced when running `plot_P_Nobs.py` and `plot_DER.py`.

## .py files
* ` generate_results.py ` Runs Monte Carlo simulations to simulate microlensing by clustered PBHs in different scenarios. One can control the value of the PBH mass, the number of PBHs per cluster, the fraction of dark matter in PBHs and the model survey. The number of microlensing events is saved in each case.
* ` plot_dgamma.py ` Plots the differential event rate for smoothly-distributed DM and chosen realisations of the Monte Carlo simulations (see Fig. 1 of the paper).
* ` plot_Nex.py ` Plots the probability distribution of the number of events against the number of events for smoothly-distributed DM and for clustered PBHs (see Fig. 2 of the paper).

# Requirements
Python 3

# License
This project is licensed under the MIT License - see the LICENSE.md file for details.

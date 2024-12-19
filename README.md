# PBH_microlensing_clustering


Code for calculating microlensing constraints for clustered primordial black hole (PBH) dark matter (DM), associated with the paper "*Effect of clustering on primordial black hole microlensing constraints*"  (see https://arxiv.org/abs/2203.04209). We use Monte Carlo simulations to estimate the probability distribution of the number of expected events, from which a 95% confidence limit on the abundance of PBHs can be estimated.


# Contents
The code, plots etc. are arranged as follows

* `data_files` Includes files required to calculate constraints for the case of smoothly-distributed DM (which has been previously considered in the literature and provides a comparison) and the differential event rate for 
smoothly-distributed DM.
* `figures` Includes figures produced when running `plot_P_Nobs.py` and `plot_DER.py`.
* ` generate_results.py ` Runs Monte Carlo simulations to simulate microlensing by clustered PBHs in different scenarios. One can control the value of the PBH mass, the number of PBHs per cluster, the fraction of dark matter in PBHs and the model survey. The number of microlensing events is saved in each case.
* ` plot_dgamma.py ` Plots the differential event rate for smoothly-distributed DM and chosen realisations of the Monte Carlo simulations (see Fig. 1 of the paper).
* ` plot_Nex.py ` Plots the probability distribution of the number of events against the number of events for smoothly-distributed DM and for clustered PBHs (see Fig. 2 of the paper).

# Requirements
Python 3

# License
This project is licensed under the MIT License - see the LICENSE.md file for details.


### UserWarning: A NumPy version >=1.16.5 and <1.23.0 is required for this version of SciPy (detected version 1.24.3


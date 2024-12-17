# PBH_microlensing_clustering


Code for calculating microlensing constraints for clustered primordial black hole (PBH) dark matter (DM), associated with the paper "*Effect of clustering on primordial black hole microlensing constraints*"  (see https://arxiv.org/abs/2203.04209). We use Monte Carlo simulations to estimate the probability distribution of the number of expected events, from which a 95% confidence limit on the abundance of PBHs can be estimated.


# Contents
The code, plots etc. are arranged as follows

* `data_files` Includes files required to calculate constraints for the case of smoothly-distributed DM (which has been previously considered in the literature and provides a comparison) and the differential event rate for 
smoothly-distributed DM.
* `Figures` Includes figures.
* ` generate_results.py ` Runs the Monte Carlo simulations. The number of microlensing events is saved in each case.
* `halomodel_optimised.py` Class describing the physical model underlying a microlensing survey, including the distance to the source stars, the model parameters used to describe the Galactic halo, the PBH mass and the number of PBHs per cluster.
* ` expected_events_discrete_clustered.py ` Includes methods for calculating the number of events from a model microlensing survey.
* ` produce_event_duration_distribution_smooth.py` Calculate the differential event rate for smoothly distributed DM, to provide a comparison to the clustered case. 
* ` compute_smooth_mean.py` Calculate the number of expected microlensing events for smoothly distributed DM, to provide a comparison to the clustered case. 
* ` plot_dgamma_v2.py ` Plots the differential event rate for smoothly-distributed DM and chosen realisations of the Monte Carlo simulations (see Fig. 1 of the paper).
* ` plot_Nex.py ` Plots the probability distribution of the number of events against the number of events for smoothly-distributed DM and for clustered PBHs (see Fig. 2 of the paper).

# Requirements
Python 2.8


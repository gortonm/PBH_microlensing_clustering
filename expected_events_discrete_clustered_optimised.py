#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 11:42:29 2021

@author: ppxmg2
"""
import numpy as np
import csv
import halomodel_optimised as hm

def load_eff():
    """
    Load efficiency function from EROS-2 for the LMC, as extracted by myself 
    from Fig. 11 of Tisserand et al. (2007)

    Returns
    -------
    eff_x : Numpy array of type Float
        Event durations, in units of days.
    eff_y : Numpy array of type Float
        Efficiency at each event duration in eff_x.

    """
    file_eff_CSV = open('./data_files/eff_func_EROS2_LMC_detailed.csv')
    eff_CSV = csv.reader(file_eff_CSV)
    list_eff_LMC = list(eff_CSV)
        
    eff_x = []
    eff_y = []
    for col in list_eff_LMC[1:]:
        eff_x.append(float(col[0]))
        # include factor of 0.9 to conservatively account for lensing by binary lenses (see caption of Fig. 9 Tisserand+ '07)
        eff_y.append(float(col[1]))
    # return efficiency function event duration values in units of years
    return np.array(eff_x) / 365.25, np.array(eff_y)


def find_tE_gamma_c(dL_values, v_values, m_pbh, n_cl, d_s=50, v_c=220):    # Find event durations and contributions to the event rate, in units of years / yr^{-1}
    setup = hm.Halomodel(m_pbh=m_pbh, n_cl=n_cl, d_s=50, v_c=220)
    return setup.einstein_radius(dL_values) / v_values, setup.event_rate(dL_values, v_values)

# calculate the number of expected events, for the discrete case
def n_ex(dL_values, v_values, m_pbh, n_cl, save=True, eff=True, EROS2=True, blendingcorrection=False, poisson=True):
    
    # calculate event durations and contributions to the total event rate
    tE_values, gamma_c_values = find_tE_gamma_c(dL_values, v_values, m_pbh, n_cl)
    
    # if set to True, boolean accounts for average blending correction from EROS-2
    if blendingcorrection:
        blend = 0.9
    else:
        blend = 1.

    # exposure values for EROS-2 observing the LMC
    if EROS2:
        n_stars = 5.49e6 # Eq. 6 Tisserand+ 2007
        obs_time = 2500 / 365.25 # convert to units of years

    # if using perfect efficiency at all event durations, return simpler expression for the number of events
    if eff == False:
        return blend * np.sum(gamma_c_values) * n_stars * obs_time

    # load efficiency function
    eff_x, eff_y = load_eff()
    
    n_obs = 0
    
    # calculate integrand at each event duration value
    for i in range(len(tE_values)):
        if poisson:
            n_obs += np.random.poisson(blend * n_stars * obs_time * gamma_c_values[i] * np.interp(tE_values[i], eff_x, eff_y, left=0, right=0) )
        else:
            n_obs += blend * n_stars * obs_time * gamma_c_values[i] * np.interp(tE_values[i], eff_x, eff_y, left=0, right=0)
                
    # return expected number of events
    return np.sum(np.array(n_obs))


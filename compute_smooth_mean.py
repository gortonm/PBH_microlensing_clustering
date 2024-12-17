#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 16:33:43 2022

@author: ppxmg2
"""

"""
Calculate the number of expected microlensing events for smoothly distributed 
DM, to provide a comparison to the clustered case.
"""

from expected_events_discrete_clustered import load_eff
import os
import numpy as np

def trapezium(x, y):
    """ Implements the trapezium rule, integrating array y with respect to 
    array x"""
    area = 0
    for i in range(len(x) - 1):
        area += (x[i+1] - x[i]) * (y[i] + 0.5*(x[i+1]-x[i])*(y[i+1]-y[i]))
    return area

m_pbh = 1000
filepath = f'{os.getcwd()}' + '/data_files/event_duration_dist_step=0.50_mpbh={:.2f}'.format(np.log10(m_pbh)) + '.txt'

t_smooth, dgamma_smooth = np.loadtxt(filepath, delimiter = ',')
t_smooth = np.array(t_smooth) / 365.25

n_stars = 5.49e6 # Eq. 6 Tisserand+ 2007
obs_time = 2500 / 365.25 # convert to units of years

eff_x, eff_y = load_eff()

mean = n_stars * obs_time * trapezium(t_smooth, dgamma_smooth * np.interp(t_smooth, eff_x, eff_y, left=0, right=0))
print(mean)



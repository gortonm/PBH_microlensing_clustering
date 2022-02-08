#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 12:31:29 2021

@author: ppxmg2
"""

import numpy as np
import halomodel_optimised as hm

def produce_values(n_cl, m_pbh, d_s, v_c, omega=84, poisson=True):
    setup = hm.Halomodel(m_pbh=m_pbh, n_cl=n_cl, d_s=d_s, v_c=v_c, omega=omega)
               
    # Choose number of samples to draw for MC simulation
    if poisson:
        n_samp = np.random.poisson(setup.n_clusters)
    else:
        n_samp = int(setup.n_clusters)

    dL_values = []
    v_values = [] 
                         
    # calculate event durations
    for j in range(n_samp):
        d_L = setup.rejection_sampler_positions()
        v = setup.sample_speeds_maxwell()
        
        dL_values.append(d_L)
        v_values.append(v)
        
    return dL_values, v_values
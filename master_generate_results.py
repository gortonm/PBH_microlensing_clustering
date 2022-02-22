# Generate samples for number of expected events

import numpy as np
import os
from expected_events_discrete_clustered import n_ex, produce_values

# Range of PBH masses to consider (in solar masses)
m_pbhs = 10**np.arange(0, 1.51, 0.5)

# Range of numbers of PBHs per cluster to consider
n_cls = 10**np.arange(3, 8.1, 1.)

# Number of realisations for each PBH mass and cluster size
n_realisations = 1000

# Set parameters for the standard halo model with the LMC at a distance of 
# 50 kpc and a circular velocity of the Sun of 220 km/s
d_s, v_c = 50, 220

def save(d_L, v, n_cl, m_pbh):
    """
    Save cluster distances and speeds (useful if I wish to produce a plot of
    the differential event rate).

    Parameters
    ----------
    d_L : Numpy array of type Float
        PBH cluster line-of-sight distances, in pc.
    v : TYPE
        PBH cluster transverse velocity, in pc / yr.
    n_cl : TYPE
        DESCRIPTION.
    m_pbh : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # Save cluster distances and speeds
    filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + str(i)
    np.savetxt(filepath + '_dL.txt', d_L)
    np.savetxt(filepath + '_v.txt', v)
    
            
# Produce a sample of cluster line of sight distances and speeds
for m_pbh in m_pbhs:
    for n_cl in n_cls:
        filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh))
        
        n_ex_EROS_efficiency = np.zeros(n_realisations)
        n_ex_perfect_efficiency = np.zeros(n_realisations)
        
        for i in range(0, n_realisations):
            np.random.seed = m_pbh * n_cl * i

            d_L, v = produce_values(n_cl, m_pbh, d_s, v_c)
            
            if n_cl >= 10**4:   # don't save cluster distances and speeds for small cluster sizes, since this uses a very large amount of storage space
                save(d_L, v, n_cl, m_pbh)
            
            # Calculate number of expected events
            n_ex_EROS_efficiency[i] = n_ex(d_L, v, m_pbh, n_cl, eff=True)    # EROS-2 efficiency curve
            n_ex_perfect_efficiency[i] = n_ex(d_L, v, m_pbh, n_cl, eff=False)    # perfect efficiency
                        
        np.savetxt(filepath + 'n_ex_EROS.txt', n_ex_EROS_efficiency)
        np.savetxt(filepath + 'n_ex_perfect.txt', n_ex_perfect_efficiency)
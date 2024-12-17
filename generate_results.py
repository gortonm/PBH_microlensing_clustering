# Generate samples for number of expected events

import numpy as np
import os

from expected_events_discrete_clustered import n_ex, produce_values, load_constraints

append = 'newRE'  # Results calculated using updated method for Einstein radius
no_save_perfect = False   # If True, don't save results assuming perfect efficiency

# Range of PBH masses to consider (in solar masses)
m_pbhs = 10**np.arange(0., 0.1, 0.5)

# Number of PBHs per cluster
n_cls = 10**np.arange(8., 6.9, -1.)

# Number of realisations for each PBH mass and cluster size
n_realisations = 10000

# Set parameters for the standard halo model with the LMC at a distance of 
# 50 kpc and a circular velocity of the Sun of 220 km/s
d_s, v_c = 50, 220


log_m_pbh_smooth, f_pbh_smooth = load_constraints(Tisserand=False)

def find_fpbh(m_pbh):
    """
    Find constraints on f_PBH for smoothly-distributed dark matter.

    Parameters
    ----------
    m_pbh : Array-like
        PBH masses, in solar masses.

    Returns
    -------
    Array-like
        Constraint on the abundance of dark matter in PBHs.

    """
    return np.interp(m_pbh, 10**np.array(log_m_pbh_smooth), f_pbh_smooth)


def setup_fpbh(m_pbh):
    """
    Produce a range of values of f_PBH around the value for smoothly-
    distributed dark matter, for input into the Monte Carlo simulations.

    Parameters
    ----------
    m_pbh : Array-like
        PBH masses, in solar masses.

    Returns
    -------
    Array-like
        Values of the fraction of dark matter in PBHs.

    """
    f_pbh_mid = round(find_fpbh(m_pbh), 2)
    return np.arange(f_pbh_mid-0.1, f_pbh_mid+0.11, 0.005)


def save(d_L, v, n_cl, m_pbh):
    """
    Save cluster distances and speeds (useful if I wish to produce a plot of
    the differential event rate).

    Parameters
    ----------
    d_L : Numpy array of type Float
        PBH cluster line-of-sight distances, in pc.
    v : Float
        PBH cluster transverse velocity, in pc / yr.
    n_cl : Float
        Number of PBHs per cluster.
    m_pbh : Float
        PBH mass, in solar masses.

    Returns
    -------
    None.

    """
    # Save cluster distances and speeds
    filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + str(i)
    np.savetxt(filepath + append + '_dL.txt', d_L)
    np.savetxt(filepath + append + '_v.txt', v)

            
# Produce a sample of cluster line of sight distances and speeds
for n_cl in n_cls:
    for m_pbh in m_pbhs:
        for f_pbh in setup_fpbh(m_pbh):
            filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh))
            
            n_ex_EROS_efficiency = np.zeros(n_realisations)
            n_ex_perfect_efficiency = np.zeros(n_realisations)

            np.random.seed(int(np.log10(m_pbh * n_cl)))
    
            for i in range(0, n_realisations):                                    
                d_L, v = produce_values(n_cl, m_pbh, d_s, v_c, f_pbh=f_pbh)
                
                if n_cl >= 10**3:   # don't save cluster distances and speeds for small cluster sizes, since this uses a very large amount of storage space
                    save(d_L, v, n_cl, m_pbh)
                
                # Calculate number of expected events
                n_ex_EROS_efficiency[i] = n_ex(d_L, v, m_pbh, n_cl)    # EROS-2 efficiency curve
                if i < 10:
                    print(n_ex_EROS_efficiency[i])
                if no_save_perfect == False:
                    n_ex_perfect_efficiency[i] = n_ex(d_L, v, m_pbh, n_cl, eff=False)    # Perfect efficiency
            
            np.savetxt(filepath + 'n_ex_EROS_2_fpbh={0:.3f}_1e4samples'.format(f_pbh) + append +'.txt', n_ex_EROS_efficiency)
            if no_save_perfect == False:
                np.savetxt(filepath + 'n_ex_perfect_fpbh={0:.3f}_1e4samples'.format(f_pbh) + append +'.txt', n_ex_perfect_efficiency)
            
np.savetxt(filepath + 'n_ex_EROS_2_fpbh={0:.3f}_1e4samples'.format(f_pbh) + append +'.txt', n_ex_EROS_efficiency)

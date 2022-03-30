# Generate samples for number of expected events

import numpy as np
import os

new_RE = False

if new_RE:
    from expected_events_discrete_clustered import n_ex, produce_values, load_constraints
else:
    from expected_events_discrete_clustered_oldRE import n_ex, produce_values, load_constraints

# Range of PBH masses to consider (in solar masses)
m_pbhs = 10**np.arange(0., 0.1, 1.)

# Number of PBHs per cluster
n_cls = 10**np.arange(8, 5.9, -1.)
#n_cls = np.array([10**3])

# Number of realisations for each PBH mass and cluster size
n_realisations = 10000

# Set parameters for the standard halo model with the LMC at a distance of 
# 50 kpc and a circular velocity of the Sun of 220 km/s
d_s, v_c = 50, 220


log_m_pbh_smooth, f_pbh_smooth = load_constraints(Tisserand=False)

def find_fpbh(m_pbh):
    return np.interp(m_pbh, 10**np.array(log_m_pbh_smooth), f_pbh_smooth)

def setup_fpbh(m_pbh):
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
    np.savetxt(filepath + 'oldRE_dL.txt', d_L)
    np.savetxt(filepath + 'oldRE_v.txt', v)
    
    
            
# Produce a sample of cluster line of sight distances and speeds
for n_cl in n_cls:
    for m_pbh in m_pbhs:
        #for f_pbh in setup_fpbh(m_pbh):
        for f_pbh in np.array([1]):
            filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh))
            
            n_ex_EROS_efficiency = np.zeros(n_realisations)
            n_ex_perfect_efficiency = np.zeros(n_realisations)
                       
            if n_cl * m_pbh > 2**32 - 1:
                np.random.seed(int(np.log(n_cl*m_pbh)))
                
            else:
                np.random.seed(int(n_cl * m_pbh))

    
            for i in range(0, n_realisations):
                
                d_L, v = produce_values(n_cl, m_pbh, d_s, v_c, f_pbh=f_pbh)
                
                
                if n_cl >= 10**3:   # don't save cluster distances and speeds for small cluster sizes, since this uses a very large amount of storage space
                    save(d_L, v, n_cl, m_pbh)
                
                # Calculate number of expected events
                n_ex_EROS_efficiency[i] = n_ex(d_L, v, m_pbh, n_cl)    # EROS-2 efficiency curve
                n_ex_perfect_efficiency[i] = n_ex(d_L, v, m_pbh, n_cl, eff=False)    # Perfect efficiency
            
            if new_RE:
                append = 'newRE'
                
            else:
                append = 'oldRE'
            np.savetxt(filepath + 'n_ex_EROS_2_fpbh={0:.3f}_1e4samples'.format(f_pbh) + append +'.txt'.format(f_pbh), n_ex_EROS_efficiency)
            np.savetxt(filepath + 'n_ex_perfect_fpbh={0:.3f}_1e4samples'.format(f_pbh) + append +'txt', n_ex_perfect_efficiency)
            

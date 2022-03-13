# Generate samples for number of expected events

import numpy as np
import os
from expected_events_discrete_clustered_areacorrection import n_ex, produce_values, event_rate_factor

def save(d_L, v, n_cl, m_pbh, i):
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
    np.savetxt(filepath + '_dL_area_factor.txt', d_L)
    np.savetxt(filepath + '_v_area_factor.txt', v)

    


# Number of PBHs per cluster
n_cl = 1000000

# Number of realisations for each PBH mass and cluster size
n_realisations = 10000

# Set parameters for the standard halo model 
d_s = 50    # source distance, in kpc
v_c = 220    # rotational speed of the Sun

# Range of PBH masses to consider (in solar masses)
#m_pbhs = np.array([1e3])
m_pbhs = 10**np.arange(0, 0.1, 1.)    # range of PBH masses to calculate for

# Range of PBH fractions to consider
f_pbhs = np.array([1])
#f_pbhs = np.arange(0.0961, 0.09651, 0.0001)    # N_cl = 1e3, M_PBH = 1e3 M_sun

#f_pbhs = np.arange(0.1, 0.2001, 0.001)    # test range (N_cl = 1e4, M_PBH = 1e3 M_sun)
#f_pbhs = np.arange(0.138, 0.1431, 0.0001)    # choice range (N_cl = 1e4, M_PBH = 1e3 M_sun)
#f_pbhs = np.arange(0.3, 0.5001, 0.001)    # test range (N_cl = 1e5, M_PBH = 1e3 M_sun)
#f_pbhs = np.arange(0.11, 0.12, 0.001) # test range (N_cl = 1e4, M_PBH = 1e2.5 M_sun))

print(len(f_pbhs))

#f_pbhs = [setup_fpbh(m_pbh) for m_pbh in m_pbhs]

#m_pbhs = 10**np.arange(2.5, -0.1, -0.5)
#f_pbhs_central_values = [0.077, 0.029, 0.014, 0.013, 0.020, 0.035, 0.076]    # values from smooth case
#f_pbhs_central_values = [0.035, 0.020, 0.013, 0.014, 0.029, 0.077]    # for m_pbhs = 10**np.arange(2.5., -0.1, -0.5)

for j, m_pbh in enumerate(m_pbhs):     # if line below is commented out, use for m_pbh in m_pbhs:
    #f_pbhs = np.arange(f_pbhs_central_values[j]-0.003, f_pbhs_central_values[j]+0.00301, 0.0001)    # set up f_PBH values, automatically chosen around smooth case
    
    """ Temporary: for saving cluster area correction on sky"""
    r_cl = 1.11198e-2 * n_cl**(5/6) * m_pbh ** (1/3)
    omega = 84 * (np.pi/180)**2

    
    for k, f_pbh in enumerate(f_pbhs):        
        
        filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh))
        
        # modify file names to record the appropriate number of samples
        if n_realisations == 20000:
            a = '_2e4samples'
        elif n_realisations == 100000:
            a = '_1e5samples'
        elif n_realisations == 200000:
            a = '_2e5samples'
        elif n_realisations == 300000:
            a = '_3e5samples'

        else:
            a = ''
        filename = filepath + 'n_ex_cutoff=200d_fpbh={0:.4f}'.format(f_pbh) + a + '_area_factor_2.txt'

        n_ex_values = np.zeros(n_realisations)
        
        seed_1 = int(n_cl*m_pbh)
        #seed_2 = 2
        #seed_3 = 3
        
        #np.random.seed(int(n_cl * m_pbh))
        
        
        """ Temporary: for saving cluster area correction on sky"""
        np.random.seed(seed_1)
        event_rate_factors = []


        for i in range(0, n_realisations):
            
            """ Temporary: for saving cluster area correction on sky"""
            d_L, v, d = produce_values(n_cl, m_pbh, d_s, v_c, f_pbh=f_pbh)
            
            for j in range(len(d_L)):
                r_cone = d_L[j] * np.sqrt(84 * (np.pi/180)**2 / np.pi)
                event_rate_factors.append(event_rate_factor(r_cl, r_cone, d[j]))
            
            # Calculate number of expected events
            n_ex_values[i] = n_ex(d_L, v, m_pbh, n_cl, d)
            
            if i < 100:
                save(d_L, v, n_cl, m_pbh, i)
            
        """ Temporary: for saving cluster area correction on sky"""
        np.savetxt(filepath + 'fcl.txt', event_rate_factors)    

        np.savetxt(filename, n_ex_values)

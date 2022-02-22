import numpy as np
import csv
import halomodel_optimised as hm

def load_eff():
    """
    Load efficiency function from EROS-2 for the LMC, extracted from Fig. 11 
    of Tisserand et al. (2007)

    Returns
    -------
    eff_x : Numpy array of type Float
        Event durations, in units of years.
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
        # Include factor of 0.9 to account for lensing by binary lenses (see caption of Fig. 9 Tisserand+ '07)
        eff_y.append(0.9*float(col[1]))
        
    return np.array(eff_x) / 365.25, np.array(eff_y)


def find_tE_gamma_c(dL_values, v_values, m_pbh, n_cl, d_s=50, v_c=220):
    """
    Find Einstein diameter crossing times (in years) and contributions to the 
    event rate (in inverse years), given values for PBH cluster distances and speeds

    Parameters
    ----------
    dL_values : Numpy array of type Float
        Array of PBH cluster line-of-sight distances, in pc.
    v_values : Numpy array of type Float
        Array of PBH cluster speeds, in pc / yr.
    m_pbh : Float
        PBH mass, in solar masses.
    n_cl : Float
        Number of PBHs per cluster.
    d_s : Float, optional
        Source distance, in kiloparsecs. The default is 50.
    v_c : Float, optional
        Circular speed in the Maxwell-Boltzmann distribution, in km/s. The default is 220.

    Returns
    -------
    Numpy array of type Float.
        Array of event durations (Einstein radius crossing times), in years.
    
    Numpy array of type Float:
        Array of contributions to the total event rate, in inverse years.

    """
    setup = hm.Halomodel(m_pbh=m_pbh, n_cl=n_cl, d_s=50, v_c=220)
    return setup.einstein_radius(dL_values) / v_values, setup.event_rate(dL_values, v_values)


def length_extended(x):
    """
    Returns length of an array, or 1 if passed a scalar of type Float.

    Parameters
    ----------
    x : Numpy array of type Float, or scalar of type Float

    Returns
    -------
    Integer
        Length of the object, 1 if a scalar of type Float.

    """
    return 1 if isinstance(x, np.float64) else len(x)


def convert_to_array(x):
    """
    Convert a scalar Float to a single-valued array, or do nothing if passed an array

    Parameters
    ----------
    x : Numpy array of type Float, or Float

    Returns
    -------
    Numpy array of type Float
        If passed scalar of type Float, return that valuue as a single-valued Numpy array of type Float.
        Otherwise, return the input array.

    """
    return np.array([x]) if isinstance(x, np.float64) else x


def produce_values(n_cl, m_pbh, d_s, v_c, omega=84, poisson=True):
    """
    Generate a sample of PBH cluster line-of-sight distances and speeds.

    Parameters
    ----------
    n_cl : Integer
        Number of PBHs per cluster.
    m_pbh : Float
        PBH mass, in solar masses.
    d_s : Float
        Microlensing source distance, in kpc.
    v_c : Float
        Circular speed of Sun in km / s.
    omega : Float, optional
        Solid angle of the viewing area for microlensing, in square degrees. The default is 84.
    poisson : Boolean, optional
        Controls whether to draw the number of PBH clusters from a Poisson distribution or use the exact value, rounded to the nearest integer. The default is True.

    Returns
    -------
    dL_values : Numpy array of type Float
        Array of PBH cluster line-of-sight distances, in pc.
    v_values : Numpy array of type Float
        Array of PBH cluster speeds, in pc / yr.

    """
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
    
    return np.array(dL_values), np.array(v_values)


def n_ex(dL_values, v_values, m_pbh, n_cl, eff=True, blendingcorrection=False, poisson=True):
    """
    Calculate the number of events in a given realisation, given values for PBH cluster distances and speeds.

    Parameters
    ----------
    dL_values : Numpy array of type Float
        Array of PBH cluster line-of-sight distances, in pc.
    v_values : Numpy array of type Float
        Array of PBH cluster speeds, in pc / yr.
    m_pbh : Float
        PBH mass, in solar masses.
    n_cl : Float
        Number of PBHs per cluster.
    eff : Boolean, optional
        If True, use the EROS-2 efficiency function. If False, assume perfect efficiency. The default is True.
    blendingcorrection : Boolean, optional
        If True, accounts for average blending correction from EROS-2 by reducing detection efficiency by 10%. The default is False.
    poisson : Boolean, optional
        Controls whether to draw the number of PBH clusters from a Poisson distribution or use the exact value, rounded to the nearest integer. The default is True.

    Returns
    -------
    Float
        Number of events in a particular realisation.

    """
    
    # calculate event durations and contributions to the total event rate
    tE_values, gamma_c_values = find_tE_gamma_c(dL_values, v_values, m_pbh, n_cl)
    
    # if set to True, boolean accounts for average blending correction from EROS-2
    if blendingcorrection:
        blend = 0.9
    else:
        blend = 1

    # Exposure values for EROS-2 observing the LMC
    n_stars = 5.49e6    # Eq. 6 Tisserand+ 2007
    obs_time = 2500 / 365.25    # convert to units of years

    # if using perfect efficiency at all event durations, return simpler expression for the number of events
    if eff == False:
        if poisson:
            return np.random.poisson(blend * np.sum(gamma_c_values) * n_stars * obs_time)
        
        else:
            return blend * np.sum(gamma_c_values) * n_stars * obs_time

    # Load efficiency function
    eff_x, eff_y = load_eff()
    
    # Number of events for a given realisation
    n_obs = 0
    
    # Ensure all inputs are Numpy arrays, even if there is only a single event in this realisation
    tE_values, gamma_c_values = convert_to_array(tE_values), convert_to_array(gamma_c_values)
    
    # Calculate integrand for number of expected events at each event duration value
    for i in range(length_extended(tE_values)):
        if poisson:
            n_obs += np.random.poisson(blend * n_stars * obs_time * gamma_c_values[i] * np.interp(tE_values[i], eff_x, eff_y, left=0, right=0) )
        
        else:
            n_obs += blend * n_stars * obs_time * gamma_c_values[i] * np.interp(tE_values[i], eff_x, eff_y, left=0, right=0)
                
    return np.sum(np.array(n_obs))


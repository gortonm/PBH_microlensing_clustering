import numpy as np
import csv
import halomodel_optimised as hm


def overlap_area(r, R, d):
    """
    Calculate the area of intersection between two overlapping circles.

    Parameters
    ----------
    r : Float
        Radius of first circle.
    R : Float
        Radius of second circle.
    d : Float
        Separation between the centres of the two circles.

    Returns
    -------
    Float
        Intersection area between the two circles.

    """
    
    if d < abs(R - r):
        return np.pi * min(r, R) **2
        
    elif d < R + r:
        if r**2 * np.arccos((d**2 + r**2 - R**2) / (2*d*r)) + R**2 * np.arccos((d**2 + R**2 - r**2) / (2 * d * R)) - 0.5*np.sqrt((-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R)) < 0:
            print('arccos1 = ',  r**2 * np.arccos((d**2 + r**2 - R**2) / (2 * d * r)))
            print('arccos2 = ', R**2 * np.arccos((d**2 + R**2 - r**2) / (2 * d * R)))
            print('sqrt term = ', 0.5*np.sqrt((-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R)))
            print('overlap area (2) = ', r**2 * np.arccos((d**2 + r**2 - R**2) / (2*d*r)) + R**2 * np.arccos((d**2 + R**2 - r**2) / (2 * d * R)) - 0.5*np.sqrt((-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R)))
            return 0
        return r**2 * np.arccos((d**2 + r**2 - R**2) / (2*d*r)) + R**2 * np.arccos((d**2 + R**2 - r**2) / (2 * d * R)) - 0.5*np.sqrt((-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R))
    
    else:
        return 0


def fractional_area_overlap(r, R, d):
    """
    Find fractional area of a circle overlapping with another circle.

    Parameters
    ----------
    r : Float
        Radius of circle to calculate fractional area for.
    R : Float
        Radius of overlapping circle.
    d : Float
        Distance between the centres of the two circles.

    Returns
    -------
    Float
        Fraction of area of circle of radius r that intersects circle of radius R.

    """
    return overlap_area(r, R, d) / (np.pi  * r**2)


def generate_d(r):
    """
    Sample distance from the centre of a circle, for uniformly distributed
    points in the circle.

    Parameters
    ----------
    r : Float
        Radius of the circle.

    Returns
    -------
    Float
        Distance from line-of-sight.

    """
    return r * np.sqrt(np.random.uniform())


def event_rate_factor(r_cl, r_cone, d):
    """
    Factor to multiply the event rate by, to account for clusters that do 
    not lie exactly along the line of sight and/or subtend a larger solid
    angle than the LMC.

    Parameters
    ----------
    d_L : Float
        PBH cluster line-of-sight distance, in pc..
    r_cl : Float
        PBH cluster radius, in pc.
    omega : Float
        Solid angle of LMC, in square radians.

    Returns
    -------
    Float
        Factor multiplying event rate, accounting for clusters that do 
        not lie exactly along the line of sight and/or subtend a larger solid
        angle than the LMC.

    """
    
    return overlap_area(r_cl, r_cone, d) / (np.pi * r_cl**2)
    


def find_tE_gamma_c(dL_values, v_values, m_pbh, n_cl, r_cl, d_values, d_s=50, v_c=220, omega=84):
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
    r_cl : Float
        PBH cluster radius, in pc.
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
    
    setup = hm.Halomodel(m_pbh=m_pbh, n_cl=n_cl, d_s=d_s, v_c=v_c, omega=omega)
    event_rate_values = np.zeros(len(dL_values))
    
    for i in range(len(dL_values)):
        
        r_cone = dL_values[i] * np.sqrt(setup.omega / np.pi)
        event_rate_values[i] = setup.event_rate(dL_values[i], v_values[i]) * event_rate_factor(r_cl=r_cl, r_cone=r_cone, d=d_values[i])

    return setup.einstein_radius(dL_values) / v_values, event_rate_values


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


def produce_values(n_cl, m_pbh, d_s, v_c, omega=84, f_pbh=1, poisson=True):
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
    
    # Create model for a given PBH cluster properties and parameters in the
    # Standard Halo Model.
    setup = hm.Halomodel(m_pbh=m_pbh, n_cl=n_cl, d_s=d_s, v_c=v_c, f_pbh=f_pbh, omega=omega)
    
    # Cluster radius (Eq. 3 of paper)
    r_cl = setup.r_cl    # checked
    
    # Radius of LMC
    r_LMC = setup.d_s * np.sqrt(setup.omega / np.pi)
    #print(r_LMC)   # checked
    
    # Ratio of volume where clusters are sampled to the volume of the microlensing cone
    volume_ratio = 1 + 3*(r_cl / r_LMC) + 3*(r_cl / r_LMC)**2
    #print(volume_ratio)    # checked
    
    # Choose number of samples to draw for MC simulation
    if poisson:
        n_samp = np.random.poisson(setup.n_clusters * volume_ratio)
    else:
        n_samp = int(setup.n_clusters * volume_ratio)

    dL_values = []
    v_values = []
    """Modified to calculate d values"""
    d_values = []
                             
    # Generate cluster distances and transverse speeds
    for j in range(n_samp):
        d_L = setup.rejection_sampler_positions()
        v = setup.sample_speeds_maxwell()
        
        r_cone = d_L * np.sqrt(setup.omega / np.pi)
        d = generate_d(r_cone + r_cl)
        
        dL_values.append(d_L)
        v_values.append(v)
        d_values.append(d)
    
    return np.array(dL_values), np.array(v_values), np.array(d_values)


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


def n_ex(dL_values, v_values, m_pbh, n_cl, d_values, poisson=True, blend=False):
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

    Returns
    -------
    Float
        Number of events in a particular realisation.

    """
    # cluster radius
    r_cl = 1.11198e-2 * n_cl**(5/6) * m_pbh ** (1/3)
    #print(r_cl)    # checked, correct
        
    # calculate event durations and contributions to the total event rate
    tE_values, gamma_c_values = find_tE_gamma_c(dL_values, v_values, m_pbh, n_cl, r_cl, d_values)
    
    # Exposure values for EROS-2 observing the LMC
    n_stars = 5.49e6    # Eq. 6 Tisserand+ 2007
    obs_time = 2500 / 365.25    # convert to units of years
    
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



def load_constraints(min_tE_200=False):
    """
    Load constraints on the PBH fraction from EROS-2.
    
    Parameters
    ----------
    min_tE_200 : Boolean. The default is False.
        If True, import constraints on the PBH fraction using a cutoff of events longer
        than 200 days.
        If False, import constraints on the PBH fraction using a cutoff of events longer
        than 100 days.

    Returns
    -------
    m_pbh : List of type Float
        PBH mass (assumed monochromatic).
    f_pbh : List of type Float
        Constraint on the dark matter fraction in PBHs at a given (monochromatic)
        PBH mass m_pbh.

    """
    if min_tE_200:
        file_constraints_CSV = open("./data_files/Constraints_min_tE_200d.csv")
    else:
        file_constraints_CSV = open("./data_files/Constraints_min_tE_100d.csv")

    constraints_CSV = csv.reader(file_constraints_CSV)
    list_constraints_CSV = list(constraints_CSV)

    m_pbh = []
    f_pbh = []
    for col in list_constraints_CSV[1:]:
        m_pbh.append(float(col[0]))
        f_pbh.append(float(col[1]))

    return m_pbh, f_pbh

# Generate samples for number of expected events

import numpy as np
import os
import warnings
import csv

# %% Set parameter values

# Fundamental constants
speed_conversion = 1.022704735e-6  # conversion factor from km/s to pc/yr
# speed of light: convert to units with [distance] = pc, [time] = yr
c = 2.99792458e5 * speed_conversion
# gravitational constant: convert to units with [distance] = pc, [time] = yr
G = 4.30091e-3 * speed_conversion ** 2

# Astrophysical parameters
d_s = 50e3  # LMC distance, in pc
v_c = 220 * speed_conversion  # Circular speed of Sun, in pc / yr
# Standard deviation for cluster transverse velocities, in pc / yr
sdev = v_c / np.sqrt(2)
# Solid angle subtended by the LMC, in square radians
omega = 84 * (np.pi / 180) ** 2
# Half-angle subtended by the LMC
theta_LMC = np.sqrt(omega/np.pi)
# Radius of region viewed in the LMC (assumed to be circular)
r_LMC = d_s * np.sqrt(omega / np.pi)

# Standard Halo Model parameters
r_0 = 8.5e3  # Distance of Sun from galactic centre, in pc
r_c = 5e3  # galactic core radius, in pc
rho_0 = 7.9e-3  # local DM density, in solar masses / pc^3
b_coord, l_coord = -32.8, 281  # LMC galactic coordinates, in degrees
# a and b parameters defined in Griest (1991)
a = (r_0 ** 2 + r_c ** 2) / d_s ** 2
b = -2 * r_0 * np.cos(np.radians(b_coord)) * np.cos(np.radians(l_coord)) / d_s

# Survey parameters
# Observing area (84 deg^2, converted to radians^2, as in EROS-2)
omega = (84 * (np.pi / 180) ** 2)
# Exposure for observing the LMC, in star yr (EROS-2)
exposure_EROS = 3.77e7
# Exposure for observing the LMC, in star yr ("toy" survey)
exposure_toy = 2.5e9


# %% User-defined parameters

# If True, use use EROS-2 efficiency function and exposure
use_eff_EROS2 = False

# If True, use 'toy' survey efficiency function and exposure
use_eff_toy = True

# if True, set PBH cluster radius to 10 pc
set_rcl_10 = True

# Cluster parameters
n_cl = 1e5  # Number of PBHs per cluster
m_pbh = 1e3  # PBH mass (in solar masses)

# Values for the fraction of dark matter in PBHs
f_pbhs = np.array([1])

# For M_PBH = 1e3 M_Sun, N_cl=1e3 (estimating 95% confidence limit on f_PBH)
# f_pbhs = np.arange(0.0945, 0.0976, 0.0001)

# Number of realisations for each PBH mass and cluster size
n_realisations = 10000     # for Figs. 1 and 2
# n_realisations = 300000     # For estimating 95% confidence limit on f_PBH

if use_eff_EROS2:
    # Load EROS-2 efficiency function
    file_eff_CSV = open('./data_files/eff_func_EROS2_LMC_detailed.csv')
    eff_CSV = csv.reader(file_eff_CSV)
    list_eff_LMC = list(eff_CSV)

    eff_x = []
    eff_y = []

    for col in list_eff_LMC[1:]:
        eff_x.append(float(col[0]))
        # Include factor of 0.9 to account for lensing by binary lenses
        # (see caption of Fig. 9 Tisserand+ '07)
        eff_y.append(0.9*float(col[1]))

# %%

# PBH cluster radius, in pc
if set_rcl_10:
    r_cl = 10

# Prefactor in the integral for the total event rate for smoothly distribute
# PBHs (see Eq. 13 of Griest 1991)
omega_0 = 2 * np.sqrt(np.pi) * rho_0 * v_c * a * \
    d_s ** (3 / 2) * np.sqrt(G) / c

# Maximum fractional distance to the LMC to include clusters from
x_max = 1

# Prefactor in the probability distribution function for cluster positions
pdf_prefactor = 1 / (
    d_s
    * (
        (
            (r_LMC ** 2 * (b ** 2 - 2 * a) - 2 * b * r_LMC * r_cl + 2 * r_cl ** 2)
            / np.sqrt(4 * a - b ** 2)
        )
        * (
            np.arctan((b + 2 * x_max) / np.sqrt(4 * a - b ** 2))
            - np.arctan(b / np.sqrt(4 * a - b ** 2))
        )
        - (
            0.5
            * r_LMC
            * (b * r_LMC - 2 * r_cl)
            * np.log((a + b * x_max + x_max ** 2) / a)
        )
        + x_max * r_LMC ** 2
    )
)

# Create path to save data files to
filepath = (f"{os.getcwd()}"
            + "/simulated_data/log10_Ncl/{0:.2f}".format(np.log10(n_cl))
            + "/log10_mpbh/{0:.2f}/".format(np.log10(m_pbh))
            )
if not os.path.exists(filepath):
    os.makedirs(filepath)
"""Note: setting a seed this way is outdated."""
np.random.seed(int(n_cl * m_pbh))

# %% Methods


def total_event_rate_integrand(x):
    """
    Compute integrand for the total event rate for smoothly distributed
    PBHs, Eq. 13 of Griest 1991.

    Parameters
    ----------
    x : Float
        Fractional line of sight distance from the observer to the LMC.

    Returns
    -------
    Float
        Integrand for the total event rate.

    """
    return np.sqrt(x * (1 - x)) / (a + b * x + x ** 2)


def total_event_rate_smooth(x, m_pbh, n_steps=10000):
    """
    Estimate total event rate for smoothly distributed PBHs, using a left
    Riemann sum to approximate Eq. 13 of Griest 1991.

    Parameters
    ----------
    x : Float
        Fractional line of sight distance from the observer to the LMC.
    m_pbh : Float
        PBH mass, in solar masses.
    n_steps : Integer, optional
        Number of terms in left Riemann sum. The default is 10000.

    Returns
    -------
    Float
        Total event rate, in yr.

    """
    x_values = np.linspace(0, x, n_steps)
    dx = x_values[1] - x_values[0]
    return (omega_0 / np.sqrt(m_pbh)) * dx * np.sum(total_event_rate_integrand(x_values))


def calculate_mass_cone_analytic(x_max=1):
    """
    Mass enclosed in microlensing cone (or in a section of the microlensing
    cone extending out to a fractional line of sight distance x_max to the
    LMC), in Solar masses.

    Parameters
    ----------
    x_max : Float
        Maximum fractional line of sight distance to include in cone (x_max=1
        includes the whole microlensing cone).
        The default is 1.

    Returns
    -------
    Float
        Total mass enclosed in microlensing cone, in solar masses.

    """
    return (f_pbh * rho_0 * a * omega * d_s ** 3) * (
        ((b ** 2 - 2 * a) / np.sqrt(4 * a - b ** 2))
        * (
            (
                np.arctan((b + 2 * x_max) / np.sqrt(4 * a - b ** 2))
                - np.arctan(b / np.sqrt(4 * a - b ** 2))
            )
        )
        - (0.5 * b * np.log((a + b * x_max + x_max ** 2) / a))
        + x_max
    )


def r_cone(d_L):
    """
    Calculate radius of microlensing cone.

    Parameters
    ----------
    d_L : Float
        Line-of-sight distance, in pc.

    Returns
    -------
    Float
        Radius of microlensing cone at d_L, in pc.

    """
    return d_L * np.sin(theta_LMC)


def generate_d(r_samp):
    """
    Sample distance from the centre of a circle, for uniformly distributed
    points in the circle.

    Parameters
    ----------
    r : Float
        Radius of the sampling region.

    Returns
    -------
    Float
        Distance of PBH cluster centre from centre of microlensing cone.
    """
    return r_samp * np.sqrt(np.random.uniform())


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
        return np.pi * min(r, R) ** 2

    elif d < R + r:
        # Occasional Runtime errors due to values in the square root being
        # calculated as negative and/or arguments of the arccos being
        # calculated with magnitude larger than one.

        # In this case, re-compute the separation of the cluster centre from
        # the line of sight and repeat.

        # When the square-root term has a very similar magnitude to the sum of
        # the two arccos() terms the computed value of area can be zero or
        # negative. This corresponds physically to the case where the
        # intersection area is tiny, but positive. In such cases, return
        # an intersection area of zero.

        try:
            area = (
                r ** 2 * np.arccos((d ** 2 + r ** 2 - R ** 2) / (2 * d * r))
                + R ** 2 * np.arccos((d ** 2 + R ** 2 - r ** 2) / (2 * d * R))
                - 0.5 * np.sqrt((-d + r + R) * (d + r - R)
                                * (d - r + R) * (d + r + R))
            )
            warnings.warn(Warning())

        except Warning:

            print("Invalid intersection area, re-calculate using a different cluster distance from centre of microlensing cone.")
            print("Cone radius R_cone [pc] = ", R)
            print("Cluster radius r_cl [pc] =", r)
            print(
                "Separation between centre of microlensing cone and cluster centre d [pc] = ", d)
            print("R+r-d [pc] =", R + r - d)

            # Re-calculate overlap area with new value of the separation between the centres of the two circles.
            d = generate_d(R + r)
            area = overlap_area(r, R, d)

        if area < 0:
            return 0

        return area

    else:
        return 0


def event_rate_factor(r_cl, r_cone, d):
    """
    Factor to multiply the event rate by, to account for clusters that do
    not lie exactly along the line of sight and/or subtend a larger solid
    angle than the LMC.

    Parameters
    ----------
    r_cl : Float
        PBH cluster radius, in pc.
    r_cone : Float
        Radius of microlensing cone, in pc.
    d : Float
        Separation of cluster centre from centre of microlensing cone, in pc.

    Returns
    -------
    Float
        Fraction of the cluster cross-sectional area that lies within the
        microlensing cone.
    """
    return overlap_area(r_cl, r_cone, d) / (np.pi * r_cl ** 2)


def einstein_radius(d_L):
    """
    Calculate Einstein radius of a lens.

    Parameters
    ----------
    d_L : Float
        Line-of-sight distance, in pc.

    Returns
    -------
    Float
        Einstein radius of a lens at line-of-sight distance d_L, in pc.

    """
    return np.sqrt(4 * G * m_pbh * d_L * (d_s - d_L) / (d_s * c ** 2))


def pdf_cluster_positions(d_L):
    """
    Probability distribution function for PBH cluster positions, from the
    extended sampling region.

    Parameters
    ----------
    d_L : Quantity
        Line-of-sight distance, in pc.

    Returns
    -------
    Float.
        Probability density of finding a cluster at line-of-sight distance
        d_L.

    """
    x = d_L / d_s
    return pdf_prefactor * (r_LMC * x + r_cl) ** 2 / (a + b * x + x ** 2)


def sample_speeds_maxwell():
    """
    Sample cluster speeds, from a Maxwellian distribution of velocities

    Parameters
    ----------


    Returns
    -------
    Float.
        Cluster transverse velocity, in pc / yr.

    """
    return np.sqrt(
        np.random.normal(loc=0, scale=sdev) ** 2
        + np.random.normal(loc=0, scale=sdev) ** 2)


def produce_values():
    """
    Generate a sample of PBH cluster line-of-sight distances and speeds.

    Parameters
    ----------

    Returns
    -------
    dL_values : Numpy array of type Float
        PBH cluster line-of-sight distances, in pc.
    v_values : Numpy array of type Float
        PBH cluster transverse velocities, in pc / yr.
    """

    # Choose number of samples to draw for MC simulation
    n_samp = np.random.poisson(n_clusters_sampling_region)

    dL_values = []
    v_values = []
    d_values = []

    # Generate cluster distances and transverse speeds
    for i in range(n_samp):
        d_L = rejection_sampler_positions()
        v = sample_speeds_maxwell()
        d = generate_d(r_cone(d_L) + r_cl)

        dL_values.append(d_L)
        v_values.append(v)
        d_values.append(d)

    return np.array(dL_values), np.array(v_values), np.array(d_values)


def save_cluster_data(dL_values, v_values, d_values, t_hat_values, event_rate_values, f_cl_values, i):
    """
    Save data on PBH clusters to a csv file.

    Parameters
    ----------
    dL_values : Numpy array of type Float
        PBH cluster line-of-sight distances, in pc.
    v_values : Numpy array of type Float
        PBH cluster transverse velocities, in pc / yr.
    d_values : Numpy array of type Float
        Separation of cluster centre from centre of microlensing cone, in pc.
    t_hat_values : Numpy array of type Float
        Array of event durations (Einstein diameter crossing times), in yr.
    event_rate_values : Numpy array of type Float
        Array of PBH cluster event rates, in yr^{-1}.
    f_cl_values : Numpy array of type Float
        Fractional area of PBH cluster within the microlensing cone.
    i : Integer
        Identifies the realisation.
    Returns
    -------
    None.

    """
    append = ""
    if set_rcl_10:
        append += "_rcl10"
    if use_eff_EROS2:
        append += "_EROS2"

    header = "D_L [yr], v [pc/yr], d [pc], Event rate [yr^{-1}], Event duration [yr], Fractional area correction"
    tosave = np.column_stack(
        (dL_values, v_values, d_values, event_rate_values, t_hat_values, f_cl_values)
    )
    np.savetxt(
        filepath + "cluster_data_" + str(i) + append + ".csv",
        tosave,
        delimiter=",",
        header=header,
        fmt="%s",
    )


def calculate_event_rate_uncorrected(d_L, v):
    """
    Calculate microlensing event rate for a single cluster, with no
    fractional area correction.

    Parameters
    ----------
    d_L : Float
        Line-of-sight distance, in pc.
    v : Float
        Cluster transverse velocity, in pc / yr.

    Returns
    -------
    Float.
        Microlensing event rate for a single cluster, in yr^{-1}, with no
        fractional area correction.
    """
    return 2 * v * n_cl * einstein_radius(d_L) / (omega * d_L ** 2)


def rejection_sampler_positions():
    """
    Sample values from the probability distribution of cluster positions,
    using rejection sampling.

    Parameters
    ----------

    Returns
    -------
    Float.
        Cluster line-of-sight distance, in pc.

    """
    inrange = False
    dL_test = np.random.uniform(0, d_s)
    while inrange == False:

        if np.random.uniform(0, pdf_max_val) <= pdf_cluster_positions(dL_test):
            inrange = True
            return dL_test
        else:
            dL_test = np.random.uniform(0, d_s)


def eff_top_hat(t_hat, t_min=400 / 365.25, t_max=15.0):
    """
    Evaluate 'top hat' efficiency function.

    Parameters
    ----------
    t_hat : Float
        Einstein diameter crossing time, in yr.
    t_min : Float, optional
        Minimum value where efficiency function is non-zero, in yr. The default is 400 / 365.25.
    t_max : Float, optional
        Maximum value where efficiency function is non-zero, in yr. The default is 15.

    Returns
    -------
    Float
        Value of 'top hat' efficiency function at event duration t_hat.
    """
    if t_min < t_hat < t_max:
        return 0.4
    else:
        return 0.0


def eff_EROS2(t_hat):
    """
    Efficiency function from EROS-2 for the LMC, extracted from Fig. 11
    of Tisserand et al. (2007).

    Parameters
    ----------
    t_hat : Float
        Einstein diameter crossing time, in years.
    Returns
    -------
    Float
        EROS-2 efficiency for events of duration t_hat.

    """
    return np.interp(365.25 * t_hat / 2, eff_x, eff_y, left=0, right=0)


def calculate_n_obs(t_hat_values, event_rate_values, use_eff_EROS2=True, use_eff_toy=False):
    """
    Calculate the number of microlensing events, given event durations, event
    event rate from each cluster and the fractional area correction for each
    cluster. Note that if both eff_EROS2 == False and eff_toy == False, the
    number of observed microlensing events is calculated assuming perfect
    efficiency with the exposure from EROS-2

    Parameters
    ----------
    t_hat_values : Numpy array of type Float
        Array of event durations (Einstein diameter crossing times), in yr.
    event_rate_values : Numpy array of type Float
        Array of PBH cluster event rates, in yr^{-1}.
    use_eff_EROS2 : Boolean
        If True, the number of observed events is calculated using the
        efficiency and exposure from EROS-2. The default is True.
    use_eff_toy : Boolean
        If True, the number of observed events is calculated using the
        efficiency and exposure from the toy survey. The default is True.

    Returns
    -------
    Float
        Number of events in a particular realisation.
    """
    n_obs = 0

    # Draw a number of events from each cluster from a Poisson distribution,
    # and find the sum to obtain the total number of events.

    for j in range(len(t_hat_values)):
        if use_eff_EROS2:
            mean = exposure_EROS * \
                event_rate_values[j] * eff_EROS2(t_hat_values[j])
        elif use_eff_toy:
            mean = exposure_toy * \
                event_rate_values[j] * eff_top_hat(t_hat_values[j])
        else:
            mean = exposure_EROS * event_rate_values[j]

        if mean.any() < 0:
            print(i)
            print("Calculated event rate [yr^{-1}] = ", event_rate_values[j])
            print("Cluster event duration [yr] = ", t_hat_values[j])
            raise Exception("Calculated mean number of events less than zero.")
        if np.isnan(mean.any()):
            print(j)
            print("Calculated event rate [yr^{-1}] = ", event_rate_values[j])
            print("Cluster event duration [yr] = ", t_hat_values[j])
            raise Exception("Calculated mean number of events is NaN.")

        else:
            n_obs += np.random.poisson(mean)

    return n_obs

# %%


# Maximum value of the cluster positions probability distribution function
pdf_max_val = pdf_cluster_positions(d_s)

# Ratio of volume where clusters are sampled to the volume of the microlensing
# cone
volume_ratio = 1 + 3 * (r_cl / r_LMC) + 3 * (r_cl / r_LMC) ** 2
if set_rcl_10:
    np.testing.assert_array_equal(round(volume_ratio, 8), 1.00666301)

for f_pbh in f_pbhs:
    print("f_PBH = ", f_pbh)

    # Mean mass of clusters in microlensing cone
    mass_cone_analytic = calculate_mass_cone_analytic()
    # Mean number of clusters in microlensing cone
    n_clusters_cone = calculate_mass_cone_analytic() / (n_cl * m_pbh)

    # Total mass in extended sampling region (truncated cone with width equal
    # to width of microlensing cone + cluster radius)
    mass_sampling_region = f_pbh * np.pi * rho_0 * a / pdf_prefactor
    print("Mass enclosed in extended sampling region [solar masses] = ",
          mass_sampling_region)

    # Mean number of clusters in sampling region
    n_clusters_sampling_region = mass_sampling_region / (n_cl * m_pbh)
    print("Mean number of clusters in extended sampling region = ",
          n_clusters_sampling_region)

    # Run to generate cluster data: line of sight distances, distances from the
    # centre of the microlensing cone, and speeds

    if "__main__" == __name__:
        append = ""
        if set_rcl_10:
            append += "_rcl10"
        if use_eff_EROS2:
            append += "_EROS2"

        filename_Nobs = (
            filepath
            + "n_ex_fpbh={0:.4f}".format(f_pbh)
            + append
            + "_nsamp={0:.1e}".format(n_realisations)
            + ".txt"
        )

        filename_Nobs_perfect = (
            filepath
            + "n_ex_perfect_fpbh={0:.4f}".format(f_pbh)
            + append
            + "_nsamp={0:.1e}".format(n_realisations)
            + ".txt"
        )

        n_obs_values = []
        n_obs_values_perfect = []
        total_event_rate_values_uncorrected = []

        for i in range(0, n_realisations):

            f_cl_values = []  # array for fractional area correction values
            dL_values, v_values, d_values = produce_values()

            for j in range(len(dL_values)):
                f_cl_values.append(event_rate_factor(
                    r_cl, r_cone(dL_values[j]), d_values[j]))

            # Calculate Einstein diameter crossing times and event rate from
            # each cluster
            t_hat_values = 2 * einstein_radius(dL_values) / v_values
            event_rate_values_uncorrected = calculate_event_rate_uncorrected(
                dL_values, v_values)

            # Save relevant data
            save_cluster_data(dL_values, v_values, d_values, t_hat_values,
                              event_rate_values_uncorrected, f_cl_values, i)

            # Include area correction factor for calculating the event rate
            event_rate_values = f_cl_values * event_rate_values_uncorrected
            n_obs_values.append(calculate_n_obs(
                t_hat_values, event_rate_values, use_eff_EROS2, use_eff_toy))
            n_obs_values_perfect.append(calculate_n_obs(
                t_hat_values, event_rate_values, use_eff_EROS2=False))

            if i % 1000 == 0:
                print('Realisation %s' % i)
                print('Number of clusters in realisation %s = ' %
                      i, len(t_hat_values))

        np.savetxt(filename_Nobs, n_obs_values)
        if use_eff_EROS2:
            np.savetxt(filename_Nobs_perfect, n_obs_values_perfect)

        print('Mean number of observed events (f_PBH = {:0.1f}) = '.format(
            f_pbh), np.mean(n_obs_values))

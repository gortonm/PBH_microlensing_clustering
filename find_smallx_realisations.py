import numpy as np
import os
from generate_results_updated_pcl_v2 import r_cone

set_rcl_10 = True
EROS2_eff = True

n_cl = 1e6  # Number of PBHs per cluster
m_pbh = 1  # PBH mass, in solar masses
speed_conversion = 1.022704735e-6  # conversion factor from km/s to pc/yr

# Astrophysical parameters
d_s = 50e3  # LMC distance, in pc

# Total number of realisations
n_realisations = 10000

# Flag samples containing clusters with a fractional line of sight distance
# smaller than this value
x_min = 0.005

realisations = []

append = ""
if set_rcl_10:
    append += "_rcl10_"
if EROS2_eff:
    append += "_EROS2_"

filepath = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + "cluster_data_update_pcl_"

def convert_to_array(x):
    """
    Convert a scalar Float to a single-valued array, or do nothing if passed an array
    Parameters
    ----------
    x : Numpy array of type Float, or Float
    Returns
    -------
    Numpy array of type Float
        If passed scalar of type Float, return that value as a single-valued Numpy array of type Float.
        Otherwise, return the input array.
    """
    return np.array([x]) if isinstance(x, np.float64) else x


for i in range(n_realisations):
        
    # Load data from csv file
    filename = filepath + str(i) + append + ".csv"    
    data = np.genfromtxt(filename, delimiter=',')
    data_transpose = np.transpose(data)
    
    # Check if array is empty
    # For sequences, (strings, lists, tuples), use the fact that empty sequences are false:
    if len(data_transpose) == 0:
        print(i)
        i += 1
        filename = f'{os.getcwd()}' + '/simulated_data_constraints/N_cl/{0:.2f}'.format(np.log10(n_cl)) + '/M_PBH/{0:.2f}/'.format(np.log10(m_pbh)) + "cluster_data_update_pcl_" + str(i) + append + ".csv"    
        data = np.genfromtxt(filename, delimiter=',')
        data_transpose = np.transpose(data)
       
    d_L = convert_to_array(data_transpose[0])
    d = convert_to_array(data_transpose[2])
    event_rate = convert_to_array(data_transpose[3])
    f_cl = convert_to_array(data_transpose[5])
    
    dL_uncorrected = d_L[d < r_cone(d_L)]
    f_cl_uncorrected = f_cl[d < r_cone(d_L)]
    event_rate_uncorrected = event_rate[d < r_cone(d_L)]
    
    if min(np.array(d_L)) / d_s < x_min:
        print('i = {:.0f}'.format(i) + ', min(x) = {:.5f}'.format(min(d_L)/50e3))
        k = np.argmin(d_L)
        print('f_cl = ', f_cl[k])
        print('D_L [pc] = ', d_L[k])
        print('Event rate / f_cl [yr^{1}] = ', event_rate[k])

        realisations.append(i)

print(realisations)
print('N_cl = {:.0e}'.format(n_cl))
print('M_PBH = {:.0e}'.format(m_pbh) + 'M_\odot')


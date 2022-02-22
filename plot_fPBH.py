# Code to plot Fig. 3 (f_PBH constraints)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import csv


# Specify the plot style
mpl.rcParams.update({"font.size": 16, "font.family": "serif"})
mpl.rcParams["xtick.major.size"] = 7
mpl.rcParams["xtick.major.width"] = 1
mpl.rcParams["xtick.minor.size"] = 3
mpl.rcParams["xtick.minor.width"] = 1
mpl.rcParams["ytick.major.size"] = 7
mpl.rcParams["ytick.major.width"] = 1
mpl.rcParams["ytick.minor.size"] = 3
mpl.rcParams["ytick.minor.width"] = 1
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["lines.linewidth"] = 1.5
mpl.rcParams["xtick.top"] = False
mpl.rcParams["ytick.right"] = False
mpl.rcParams["font.family"] = "serif"
mpl.rc("text", usetex=True)

mpl.rcParams["legend.edgecolor"] = "inherit"


def find_percentiles(x, minval=0.16, maxval=0.84):
    """
    Find percentile values from a data set.

    Parameters
    ----------
    x : Numpy array of type Float
        Input data to find upper and lower percentiles of.
    minval : Float, optional
        Lowest percentile to calculate. The default is 0.16.
    maxval : Float, optional
        Highest percentile to calculate. The default is 0.84.

    Returns
    -------
    Float
        Percentile of x corresponding to minval.
    Float
        Percentile of x corresponding to maxval.

    """
    x = np.sort(x)
    return x[int(len(x) * minval) - 1], x[int(len(x) * maxval) - 1]


def error_bars(x):
    """
    Find upper and lower 68% and 95% percentiles from a data set.

    Parameters
    ----------
    x : Numpy array of type Float.
        Input data to find 68% and 95% upper and lower percentiles of.

    Returns
    -------
    List of type Float
        16% and 84% percentile in x.
    List of type Float
        2.5% and 97.5% percentile in x.

    """
    x_68_lower, x_68_upper = find_percentiles(x)
    x_95_lower, x_95_upper = find_percentiles(x, minval=0.025, maxval=0.975)
    
    return [[np.median(x) - x_68_lower], [x_68_upper - np.median(x)]], [[np.median(x) - x_95_lower], [x_95_upper - np.median(x)]]


def load_constraints(Tisserand=True):
    """
    Load constraints on the PBH fraction from EROS-2, from Fig. 11 of Tisserand et al. (2007) (green curve)

    Returns
    -------
    m_pbh : List of type Float
        PBH mass (assumed monochromatic).
    f_pbh : List of type Float
        Constraint on the dark matter fraction in PBHs at a given (monochromatic)
        PBH mass m_pbh.

    """
    if Tisserand:
        name = "tisserand"
    else:
        name = "green_reproduction"

    file_constraints_CSV = open("./data_files/eros2_constraints_" + name + ".csv")
    constraints_CSV = csv.reader(file_constraints_CSV)
    list_constraints_CSV = list(constraints_CSV)

    m_pbh = []
    f_pbh = []
    for col in list_constraints_CSV[1:]:
        m_pbh.append(float(col[0]))
        f_pbh.append(float(col[1]))

    return m_pbh, f_pbh


# Range of PBH masses to consider
m_pbhs = 10 ** np.arange(0.0, 1.51, 0.5)

# Range of numbers of PBHs per cluster to consider
n_cls = 10 ** np.arange(3, 8.1, 1.0)

# Number of realisations for each PBH mass and cluster size
n_realisations = 1000

# Colours to plot for different cluster sizes
colours = ["darkturquoise", "b", "m", "darkorange", "r", "saddlebrown"]


plt.figure(figsize=(8, 4))

# Plot constraint for smoothly-distributed PBHs from EROS-2 (Fig. 15 of Tisserand et al. 2007)
m_pbh_smooth, f_pbh_smooth = load_constraints()
plt.plot(m_pbh_smooth, f_pbh_smooth, color="lime")

# Plot reproduced constraint for smoothly-distributed PBHs from EROS-2 (Fig. 1 of Green 2016)
m_pbh_smooth, f_pbh_smooth = load_constraints(Tisserand=False)
plt.plot(m_pbh_smooth, f_pbh_smooth, linestyle="dotted", color="k")

for j, n_cl in enumerate(n_cls):

    # Displace positions of different cluster sizes on the x-axis
    dx = 0.035 * (-len(n_cls) / 2 + j * 1.2)

    for m_pbh in m_pbhs:

        # Load file for number of events in a particular realisation
        filepath = (
            f"{os.getcwd()}"
            + "/simulated_data_constraints/N_cl/{0:.2f}".format(np.log10(n_cl))
            + "/M_PBH/{0:.2f}/".format(np.log10(m_pbh))
        )
        n_ex_EROS_efficiency = np.loadtxt(filepath + "n_ex_EROS.txt")

        f_pbhs = np.ones(len(n_ex_EROS_efficiency))
        for i in range(len(n_ex_EROS_efficiency)):
            if n_ex_EROS_efficiency[i] > 3:
                f_pbhs[i] = 3 / n_ex_EROS_efficiency[i]

        # Calculate arrays of 68% and 95% upper and lower percentiles
        err_68, err_95 = error_bars(f_pbhs)

        # Plot percentile range bars
        plt.errorbar(
            np.log10(m_pbh) + dx,
            np.median(f_pbhs),
            yerr=err_68,
            marker="x",
            color=colours[j],
            elinewidth=1,
            capsize=2,
        )
        plt.errorbar(
            np.log10(m_pbh) + dx,
            np.median(f_pbhs),
            yerr=err_95,
            marker="x",
            color=colours[j],
            elinewidth=0.5,
            capsize=1,
        )

    # Plot median value of f_PBH
    point = plt.plot(
        np.log10(m_pbh) + dx,
        np.median(f_pbhs),
        color=colours[j],
        marker="x",
        linestyle="none",
        label="$10^{:.0f}$".format(np.log10(n_cl)),
    )

plt.legend(title="$N_\mathrm{cl}$")

plt.xlim(-0.4, 2.1)
plt.ylim(0.0, 1.01)
plt.xlabel("$\log_{10}(M_\mathrm{PBH}/M_\odot)$")
plt.ylabel("$f_\mathrm{PBH}$")
plt.tight_layout()
plt.savefig("f_PBH_MPBH.pdf")

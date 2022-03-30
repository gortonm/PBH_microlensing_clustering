# Class for a particular set-up

import numpy as np


class Halomodel:
    def __init__(
        self,
        m_pbh=1,
        n_cl=10,
        d_s=50,
        v_c=220,
        f_pbh=1.0,
        b_coord=-32.8,
        l_coord=281,
        omega=84,
        rho_0=7.9e-3,
        r_c=5e3,
        r_0=8.5e3,
    ):
        """

        Parameters
        ----------
        m_pbh : Float, optional
            PBH mass, in solar masses. The default is 1.
        n_cl : Float, optional
            Number of PBHs per cluster. The default is 10.
        d_s : Float, optional
            Source distance, in kiloparsecs. The default is 50.
        v_c : Float, optional
            Circular speed in the Maxwell-Boltzmann distribution, in km/s.
            The default is 220.
        f_pbh : Float, optional
            Fraction of dark matter in primordial black holes. The default is 1.
        b_coord : Float, optional
            Large Magellanic Cloud latitude (in galactic coordinates), in degrees. The default is -32.8.
        l_coord: Float, optional
            Large Magellanic Cloud longitude (in galactic coordinates), in degrees. The default is 281.
        omega : Float, optional
            Solid angle of the viewing area for microlensing, in square degrees. The default is 84.
        rho_0 : Float, optional
            Local dark matter density, in solar masses / pc^3. The default is 7.9e-3.
        r_c : Float, optional
            Core radius in the standard halo model, in pc. The default is 5e3.
        r_0 : Float, optional
            Distance of the Sun from the galactic centre, in pc. The default is 8.5e3.

        Returns
        -------
        None.

        """
        # conversion factor from km/s to pc/yr
        speed_conversion = 1.022704735e-6
        
        self.c, self.G = 2.99792458e5 * speed_conversion, 4.30091e-3 * speed_conversion**2     # convert to units with [distance] = pc, [time] = yr


        self.m_pbh = m_pbh  # mass of each PBH
        self.f_pbh = f_pbh  # Fraction of DM in PBHs
        self.omega = (
            omega * (np.pi / 180) ** 2
        )  # 84 deg^2 observing area, converted to radians^2, for EROS-2 (Tisserand+ 2007 Abstract)
        self.n_cl = n_cl  # number of PBHs in each cluster
        self.v_c = v_c * speed_conversion  # circular velocity of the Sun

        # Standard deviation for cluster transverse velocities
        self.sdev = self.v_c / np.sqrt(2)

        # Set astrophysical parameters
        self.rho_0 = rho_0  # dark matter density at the Sun, in units M_sun pc^{-3}
        self.r_c = r_c  # galaxy core radius, in pc
        self.r_0 = r_0  # Sun's distance from galactic centre, in pc
        self.b_coord, self.l_coord = (
            b_coord,
            l_coord,
        )  # galactic coordinates, in degrees
        self.d_s = d_s * 1e3  # distance of source, in pc (using Griest 1991 value)

        # Parameters appearing in probability distribution function for cluster positions
        self.a = (
            -2
            * self.r_0
            * np.cos(np.radians(self.b_coord))
            * np.cos(np.radians(self.l_coord))
        )
        self.b = self.r_0**2 + self.r_c**2
        self.k = (
            0.5
            * np.sqrt(4 * self.b - self.a**2)
            * (
                (
                    np.arctan(
                        (self.a + 2 * self.d_s) / np.sqrt(4 * self.b - self.a**2)
                    )
                    - np.arctan(self.a / np.sqrt(4 * self.b - self.a**2))
                )
            )
        ) ** (-1)

        # Mean number of clusters in a conical volume, for the standard halo model and a monochromatic PBH mass function.
        self.n_clusters = (
            self.rho_0 * self.b * self.f_pbh * self.omega / (self.n_cl * self.m_pbh)
            ) * (
            (
                ((self.a ** 2 - 2 * self.b) / np.sqrt(4 * self.b - self.a ** 2))
                * (
                    (
                        np.arctan(
                            (self.a + 2 * self.d_s) / np.sqrt(4 * self.b - self.a ** 2)
                        )
                        - np.arctan(self.a / np.sqrt(4 * self.b - self.a ** 2))
                    )
                )
            )
            - 0.5
            * (self.a
            * np.log(((self.d_s * (self.a + self.d_s)) + self.b) / self.b))
            + self.d_s
        )
                
        # Normalisation factor for cluster positions probability distribution function
        self.k_cl = (
            self.omega
            * self.rho_0
            * self.b
            * self.f_pbh
            / (self.n_clusters * self.n_cl * self.m_pbh)
        )

        # Maximum value of the cluster positions probability distribution function
        self.pdf_max_val = np.max(
            [
                self.pdf_cluster_positions(-2 * self.b / self.a),
                self.pdf_cluster_positions(self.d_s),
            ]
        )

        # cluster radius, in pc
        self.r_cl = 1.11198e-2 * self.n_cl**(5/6) * self.m_pbh ** (1/3)
        
        # limiting value of cluster distance below which cluster subtends a larger solid angle
        # than the cluster itself
        self.dL_lim = self.r_cl * np.sqrt(np.pi / self.omega)
                

    def einstein_radius(self, d_L):
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
        return 4.371625683e-7 * np.sqrt(
            self.m_pbh * (self.d_s - d_L) * (d_L / self.d_s)
        )

    def event_rate(self, d_L, v_c):
        """
        Calculate microlensing event rate for a single cluster.

        Parameters
        ----------
        d_L : Float
            Line-of-sight distance, in pc.
        v_c : Quantity
            Cluster transverse velocity, in pc / yr.

        Returns
        -------
        Float.
            Microlensing event rate for a single cluster, in yr^{-1}.

        """
        return 2 * v_c * self.n_cl * self.einstein_radius(d_L) / (self.omega * d_L**2)
        """
        if d_L > self.dL_lim:
            return 2 * v_c * self.n_cl * self.einstein_radius(d_L) / (self.omega * d_L**2) 
        else:
            factor = self.omega / (np.pi * self.r_cl**2 / d_L**2)
            return factor * 2 * v_c * self.n_cl * self.einstein_radius(d_L) / (self.omega * d_L**2)
        """
    def pdf_cluster_positions(self, d_L):
        """
        Probability distribution function for PBH cluster positions.

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

        return (
            (
                self.omega
                * self.rho_0
                * self.b
                * self.f_pbh
                / (self.n_clusters * self.n_cl * self.m_pbh)
            )
            * d_L**2
            / (d_L**2 + (self.a * d_L) + self.b)
        )

    def sample_speeds_maxwell(self):
        """
        Sample cluster speeds, from a Maxwellian distribution of velocities

        Parameters
        ----------


        Returns
        -------
        Float.
            Cluster transverse velocity, in units of pc / yr.

        """
        return np.sqrt(
            np.random.normal(loc=0, scale=self.sdev) ** 2
            + np.random.normal(loc=0, scale=self.sdev) ** 2
        )

    def rejection_sampler_positions(self):
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
        x_test = np.random.uniform(0, self.d_s)
        while inrange == False:

            if np.random.uniform(0, self.pdf_max_val) <= self.pdf_cluster_positions(
                x_test
            ):
                inrange = True
                return x_test
            else:
                x_test = np.random.uniform(0, self.d_s)

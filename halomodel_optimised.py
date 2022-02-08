#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 15:51:17 2021

@author: ppxmg2
"""

import numpy as np

class Halomodel:
    def __init__(self, m_pbh=1, n_cl=10, d_s=50, v_c=220, f_pbh=1.0, b_coord = -32.8, l_coord = 281, omega=84, rho_0=7.9e-3, a_0=5e3, r_0=8.5e3, n_stars=1.*10**7, Renault98=False, test=False):
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
        n_stars : Float, optional
            Number of stars observed.
        test : Boolean, optional
            If True, use parameters of the pdf that allow for simple analytic
            results to test the code.

        Returns
        -------
        None.

        """
        # conversion factor from km/s to pc/yr
        #speed_conversion = 1.022704735e-6
        speed_conversion = 1.0227e-6

        self.m_pbh = m_pbh   # mass of each PBH
        self.f_pbh = f_pbh  # Fraction of DM in PBHs        
        self.omega = omega * (np.pi / 180) ** 2  # 84 deg^2 observing area, converted to radians^2, for EROS-2 (Tisserand+ 2007 Abstract)

        self.n_cl = n_cl  # number of PBHs in each cluster

        self.v_c = v_c * speed_conversion # Circular velocity

        # standard deviation for speeds
        self.sdev =self.v_c / np.sqrt(2)
        
        """ Set astrophysical parameters"""
        self.rho_0 = rho_0   # dark matter density at the Sun, in units M_sun pc^{-3}
        self.a_0 = a_0   # galaxy core radius, in pc
        self.r_0 = r_0   # Sun's distance from galactic centre, in pc
        self.b_coord, self.l_coord = b_coord, l_coord  # galactic coordinates, in degrees
        self.d_s = (
            d_s * 1e3 
        )  # distance of source, in pc (using Griest 1991 value)
        self.a = (
            -2
            * self.r_0
            * np.cos(np.radians(self.b_coord))
            * np.cos(np.radians(self.l_coord))
        )
        self.b = self.r_0 ** 2 + self.a_0 ** 2
        self.k = (
            0.5
            * np.sqrt(4 * self.b - self.a ** 2)
            * (
                (
                    np.arctan(
                        (self.a + 2 * self.d_s) / np.sqrt(4 * self.b - self.a ** 2)
                    )
                    - np.arctan(self.a / np.sqrt(4 * self.b - self.a ** 2))
                )
            )
        ) ** (-1)
        
        # (Mean) number of clusters in a conical volume, for the standard halo model and a monochromatic PBH mass function.
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
                
        # normalisation factor for cluster positions PDF
        self.k_cl = (
                self.omega
                * self.rho_0
                * self.b
                * self.f_pbh
                / (self.n_clusters * self.n_cl * self.m_pbh)
            )

                
        # max value of the cluster line of sight distance PDF
        self.pdf_max_val = np.max(
            [self.pdf_cluster_positions(-2 * self.b / self.a) , self.pdf_cluster_positions(self.d_s) ]
        )

    def einstein_radius(self, d_L):
        """
        Calculate Einstein radius of a lens.

        Parameters
        ----------
        d_L : Quantity
            Line-of-sight distance, in pc.

        Returns
        -------
        Float
            Einstein radius of a lens at line-of-sight distance d_L, in pc.

        """
        #return 4.371625683e-7 * np.sqrt(self.m_pbh * (self.d_s - d_L) * (d_L / self.d_s))
        return 4.37e-7 * np.sqrt(self.m_pbh * (self.d_s - d_L) * (d_L / self.d_s))


    def event_rate(self, d_L, v_c, u_T=1):
        """
        Calculate microlensing event rate for a single cluster.

        Parameters
        ----------
        d_L : Float
            Line-of-sight distance, in pc.
        v_c : Quantity
            Cluster velocity, in pc / yr.
        u_T : Float, optional
            Threshold impact parameter, in units of the Einstein radius.
            The default is 1.

        Returns
        -------
        Float.
            Microlensing event rate for a single cluster.

        """
        return  (2 * v_c * self.n_cl * self.einstein_radius(d_L) * u_T / (self.omega * d_L ** 2))
            
            
    def pdf_cluster_positions(self, d_L):
        """
        Probability distribution function for PBH cluster positions.

        Parameters
        ----------
        d_L : Quantity
            Line-of-sight distance, with dimensions of length.

        Returns
        -------
        Float.
            Probability density of finding a cluster at line-of-sight distance
            d_L.

        """
        
        """
        if d_L >= self.d_s or d_L < 0:
            return 1e-16
        else:
           return (self.omega
           * self.rho_0
           * self.b
           * self.f_pbh
           / (self.n_clusters * self.n_cl * self.m_pbh)
           ) * d_L ** 2 / (d_L ** 2 + (self.a * d_L) + self.b)
        """
        return (self.omega
           * self.rho_0
           * self.b
           * self.f_pbh
           / (self.n_clusters * self.n_cl * self.m_pbh)
           ) * d_L ** 2 / (d_L ** 2 + (self.a * d_L) + self.b)
        
    
    def optical_depth(self, d_L):
        """
        Calculate the optical depth due to a single cluster.

        Parameters
        ----------
        d_L : Quantity
           Line of sight distance, with dimensions of length.

        Returns
        -------
        Float.
            Contribution to the total optical depth from a single cluster.

        """
        return self.n_cl * np.pi * self.einstein_radius(d_L)**2 / (self.omega * d_L**2)


    def sample_speeds_maxwell(self):
        """
        Sample cluster speeds, from a Maxwellian distribution of velocities

        Parameters
        ----------


        Returns
        -------
        Float.
            Sample cluster speed, in units of pc / yr.

        """
        return np.sqrt(np.random.normal(loc=0, scale=self.sdev)**2 + np.random.normal(loc=0, scale=self.sdev)**2)
    
    
    def rejection_sampler_positions(self):
        """
        Sample values from a probability distribution using rejection sampling.
        
        
        Parameters
        ----------
    
        Returns
        -------
        Float.
            Sample from the probability distribution.
    
        """
        inrange = False
        x_test = np.random.uniform(0, self.d_s)
        while inrange == False:
    
            if np.random.uniform(0, self.pdf_max_val) <= self.pdf_cluster_positions(x_test):
                inrange = True
                return x_test
            else:
                x_test = np.random.uniform(0, self.d_s)
            

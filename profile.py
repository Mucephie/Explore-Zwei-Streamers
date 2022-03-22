
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import units as un
from astropy import constants as cnst

## additional imports
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) 

import gala.potential as gp
from gala.units import galactic

from galpy.potential import MWPotential2014
from galpy.potential import evaluateDensities
from galpy.util import conversion




'''
# def profileGala(m, r, r_s):
#     # convert to solar masses
#     rho = - cnst.G * m * np.log10(1 + r / r_s) / r
#     # {'x': x, 'y': y, 'z': z}, {'m': m, 'r_s': r_s, 'a': a, 'b': b, 'c': c, 'G': G})
#     return rho


# def profileMcMillan(rho_s, r, r_s, gam = 1):
#     # gamma is 1 for NFW according to paper
#     rho = rho_s / ((r / r_s)**gam * (1 + (r / r_s))**(3 - gam) )
#     return rho

# def velocity_nfw(v_h, r, r_s):
#     # lin n' li
#     v_nfw = v_h * np.sqrt((r_s / r) * np.log((r + r_s) / r_s) - (r_s / (r + r_s)))
#     return v_nfw

# def profileEinasto(rho_s, n, r, r_s):
#     # Jiao
#     rho = rho_s * np.exp(- (r / r_s)**(1 / n))
#     return rho

# Sofue

# Bonaca 8.1kpc 220kms --> Jerabkova 12kpc 1.07E11 2.565kpc 0.24E11
'''




class profile:
    def __init__(self, r_s):
        self.r_s = r_s




class GalaProfile():
    def __init__(self, r_s, m):
        self.r_s = r_s
        self.m = m
    def density(self, r):
        # convert to solar masses
        return - cnst.G * self.m * np.log10(1 + r / self.r_s) / r
    
class McMillanProfile():
    def __init__(self, r_s, rho_s, gamma = 1, v_h = None):
        self.r_s = r_s
        self.rho_s = rho_s
        self.gamma = gamma
        self.v_h = v_h
    def density(self, r):
        return self.rho_s / ((r / self.r_s)**self.gamma * (1 + (r / self.r_s))**(3 - self.gamma))
    def velocity(self, r, v_h = None):
        if v_h == None:
            if self.v_h == None:
                raise Exception('Error, No scale velocity provided.')
            else:
                v_h = self.v_h
        else:
            self.v_h = v_h
        return v_h * np.sqrt(((self.r_s / r) * np.log(((r + self.r_s) / self.r_s).value) - (self.r_s / (r + self.r_s))).value)

class EinastoProfile():
    def __init__(r_s, self, rho_s, n):
        self.r_s = r_s
        self.rho_s = rho_s
        self.n = n
    def density(self, r):
        return self.rho_s * np.exp(- (r / self.r_s).value**(1 / self.n))

class GalaProfileWrap():
    def __init__(r_s, self, m):
        self.r_s = r_s
        self.m = m
        self.potential = gp.NFWPotential(m = self.m, r_s = self.r_s, units = galactic)
        self.halpy = self.potential.to_galpy_potential()
    def density(self, r, use_gala = False):
        if use_gala:
            return self.potential.density(r)
        else:
            dens = []
            for i in range(len(r)):
                dens.append(self.halpy.dens(r, 0).value)
            return dens
    
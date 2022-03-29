
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

import sympy
from sympy import symbols

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


rho_crit = 3 * cosmo.H(0)**2 / (8 * np.pi * cnst.G)
# print(rho_crit.to(un.kg / un.m**3))

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
    
class LumpyProfile():
    def __init__(self, r_s, rho_s, gamma = 1):
        self.r_s = r_s
        self.rho_s = rho_s
        self.gamma = gamma
        # number density and subhalo mass
        # subhalo distribution
        # halo mass vs. total subhalo mass


    def density(self, r):
        return self.rho_s / ((r / self.r_s)**self.gamma * (1 + (r / self.r_s))**(3 - self.gamma))
    




class NFWProfile():
    def __init__(self, r_s, rho_s):
        self.r_s = r_s
        self.rho_s = rho_s # properly define this, should I use rho_0?
       
        self.M = None
        self.mass_profile = None

        self.potential_profile = None

        # self.c = 10 # 10 to 15 for mw via wikipedia, concentration parameter
        # self.R_vir = None # R_vir = c * r_s 
        # self.r_max = None

        self.r_200 = None
        self.m_200 = None

    def density(self, r):
        return self.rho_s / ((r / self.r_s)**1 * (1 + (r / self.r_s))**(3 - 1))

    def mass(self, r = None, set = False):
        if r == None:
            if self.r_200 != None:
                r = self.r_200
            else:
                self.set_200()
                r = self.r_200
        if set:
            self.mass_profile = 4 * np.pi * self.rho_s * self.r_s**3 * (np.log((self.r_s + r) / self.r_s)  + (self.r_s / (self.r_s + r) - 1))
        else:
            return 4 * np.pi * self.rho_s * self.r_s**3 * (np.log((self.r_s + r)/self.r_s)  + (self.r_s/(self.r_s + r) - 1))

    def r_when(self, rho): # TODO
        r= symbols('r', real = True)
        expr = (self.rho_s / rho).value / ((r) * (1 + (r ))**(2)) # ((sympy.sqrt(r)*(1 + r**2)) / sympy.sqrt((self.rho_s / rho).value))
        print('soving...')
        rr_s = sympy.solve(sympy.Eq(expr, 1), r)
        # print('ratio of r/r_s: ', rr_s) # rr_s = r / self.r_s
        return rr_s * self.r_s
        
    def set_200(self):
        rho_crit = (3 * cosmo.H(0)**2 / (8 * np.pi * cnst.G)).to(un.Msun / un.pc**3)
        self.r_200 = self.r_when(200 * rho_crit).to(un.kpc)
        self.m_200 = (100 * self.r_200**3 * cosmo.H(0)**2 / cnst.G).to(un.Msun)

    def potential(self, r, set = False):
        if set:
            self.potential_profile = ((4 * np.pi * cnst.G * self.rho_s * self.r_s**3) / r) * np.log(1 + (r / self.r_s))
        else:
            return ((4 * np.pi * cnst.G * self.rho_s * self.r_s**3) / r) * np.log(1 + (r / self.r_s))
    
    
    ## r_200 is when density is 200*rho_crit
    ## m_200 is 100 * r_200**3 * H**2 / cnst.G

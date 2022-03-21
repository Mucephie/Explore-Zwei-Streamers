import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import units as un
from astropy import constants as cnst


def profileGala(m, r, r_s):
    # convert to solar masses
    rho = - cnst.G * m * np.log10(1 + r / r_s) / r
    # {'x': x, 'y': y, 'z': z}, {'m': m, 'r_s': r_s, 'a': a, 'b': b, 'c': c, 'G': G})
    return rho


def profileMcMillan(rho_s, r, r_s, gam = 1):
    # gamma is 1 for NFW according to paper
    rho = rho_s / ((r / r_s)**gam * (1 + (r / r_s))**(3 - gam) )
    return rho

def velocity_nfw(v_h, r, r_s):
    # lin n' li
    v_nfw = v_h * np.sqrt((r_s / r) * np.log((r + r_s) / r_s) - (r_s / (r + r_s)))
    return v_nfw

def profileEinasto(rho_s, n, r, r_s):
    # Jiao
    rho = rho_s * np.exp(- (r / r_s)**(1 / n))
    return rho

# Sofue

# Bonaca 8.1kpc 220kms --> Jerabkova 12kpc 1.07E11 2.565kpc 0.24E11
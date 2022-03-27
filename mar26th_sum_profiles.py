## imports ##
import profile as pf

import numpy as np
import matplotlib.pyplot as plt

from astropy import units as un
from astropy import constants as cnst

import matplotlib.collections as collections
from matplotlib.patches import Polygon

import gala.potential as gp
from gala.units import galactic

from galpy.potential import MWPotential2014
from galpy.potential import evaluateDensities
from galpy.util import conversion


## Params ##
r_s   = 8.8  * un.kpc # 8.8 , 8.1 Pal5/Schonrich2010
v_h   = 220  * un.km / un.s # 220 Schonrich2010
m     = 6E11 * un.Msun
rho_s = 0.11    * un.Msun / (un.pc)**3 # 0.11
r     = np.logspace(-3, 1.5, 256) * un.kpc

## functions ##
def mkProfPlot(profile, r, label, a = 14, b = 21, shade = False):
    if shade:
        # Make the shaded region
        ix = np.linspace(a, b) * un.kpc
        iy = profile.density(ix)
        verts = [(a, 0), *zip(ix.value, iy.value), (b, 0)]
        poly = Polygon(verts, facecolor='0.5', edgecolor='0.5')
        ax.add_patch(poly)
        verts = [(a, 10E3), *zip(ix.value, iy.value), (b, 10E3)]
        poly = Polygon(verts, facecolor='0.9', edgecolor='0.5')
        ax.add_patch(poly)
    ax.semilogy(r, profile.density(r), label = label) # semilogx

## Milky way potential 2014
densmw = [] # evaluateDensities(MWPotential2014,1.,0.)*conversion.dens_in_msolpc3(220.,8.)
densmw0, densmw1, densmw2= [], [], []
for i in range(len(r)):
    x = r[i]

    # densmw.append(evaluateDensities(MWPotential2014,r, 0.)*conversion.dens_in_msolpc3(vel_r,r))
    densmw0.append(MWPotential2014[0].dens(x, 0.))
    densmw1.append(MWPotential2014[1].dens(x, 0.))
    densmw2.append(MWPotential2014[2].dens(x, 0.))
    densmw.append(evaluateDensities(MWPotential2014,x, 0.))

# print(MWPotential2014[1].dens(15, 0.))



## initialize June's profiles 2022 ##
nfwp  = pf.NFWProfile(r_s = r_s, rho_s = rho_s)
# nfwp  = pf.McMillanProfile(r_s = r_s, rho_s = rho_s, v_h = v_h)
nfwp_v2  = pf.McMillanProfile(r_s = 8.1 * un.kpc, rho_s = rho_s, v_h = v_h)
nfwp_v3  = pf.McMillanProfile(r_s = 15 * un.kpc, rho_s = 0.011 * un.Msun / (un.pc)**3 , v_h = v_h)

# print(nfwp.density(15 * un.kpc))
# print(nfwp_v2.density(15 * un.kpc))
# print(nfwp_v3.density(15 * un.kpc))


# print('r_when: ', nfwp.r_when(0.1 * un.Msun / (un.pc)**3))









## Graph it ##
fig, ax = plt.subplots()
ax.axvline(15, color='black', lw=0.5)
ax.axvline(25, color='black', lw=0.5)
# ax.axvline(0.946008635110769)


# profiles plot #
mkProfPlot(nfwp, r, 'NFW_8.8kpc', shade = True)
# mkProfPlot(nfwp_v2, r, 'NFW_8.1kpc')
# mkProfPlot(nfwp_v3, r, 'NFW_15kpc')
ax.semilogy(r, densmw, label = 'mw')

# anatomy of a plot #
plt.legend()
ax.set_title('DM Halo Mass Density')
ax.set_xlabel(r'radius ($kpc$)')
ax.set_ylabel(r'mass density ($M_{sol}/{\rm pc}^{3}$)')
ax.set_ylim((0, 10E3))
ax.set_xlim((0, 10**1.5))
plt.grid()
plt.show()
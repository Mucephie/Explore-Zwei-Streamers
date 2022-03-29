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

from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) 

## Params ##
r_s   = 8.9  * un.kpc # 8.8 , 8.1 Pal5/Schonrich2010
v_h   = 220  * un.km / un.s # 220 Schonrich2010
m     = 6E11 * un.Msun
rho_s = 0.01 * un.Msun / (un.pc)**3 # 0.11
r     = np.logspace(-3, 1.5, 256) * un.kpc
# 6.0×104ρcrit
rho_crit = 3 * cosmo.H(0)**2 / (8 * np.pi * cnst.G)
print((6*10**4*rho_crit).to(un.Msun / un.pc**3))
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


## initialize June's profiles 2022 ##
nfwp  = pf.NFWProfile(r_s = r_s, rho_s = rho_s)
nfwp.mass(r, True)
nfwp.set_200()
print('mass: ', (nfwp.m_200 / 10**12).value[0], ' x 10^12 Msun, radius: ',  nfwp.r_200.value[0], ' kpc')
nfwp.potential(r, True)


## Milky way potential 2014
densmw = [] # evaluateDensities(MWPotential2014,1.,0.)*conversion.dens_in_msolpc3(220.,8.)
densmw0, densmw1, densmw2, densNFW, denssum = [], [], [], [], []
for i in range(len(r)):
    x = r[i]

    # densmw.append(evaluateDensities(MWPotential2014,r, 0.)*conversion.dens_in_msolpc3(vel_r,r))
    densmw0.append(MWPotential2014[0].dens(x, 0.))
    densmw1.append(MWPotential2014[1].dens(x, 0.))
    densmw2.append(MWPotential2014[2].dens(x, 0.))
    densNFW.append(nfwp.density(r[i]))
    densmw.append(evaluateDensities(MWPotential2014,x, 0.))

# print(densNFW)
for i in range(len(r)):
    denssum.append(densNFW[i].value + densmw0[i] + densmw1[i])
# print(MWPotential2014[1].dens(15, 0.))






# Graph density ##
fig, ax = plt.subplots()
ax.axvline(15, color='black', lw=0.5)
ax.axvline(25, color='black', lw=0.5)
# ax.axvline(0.946008635110769)


# profiles plot #
mkProfPlot(nfwp, r, 'NFW_8.8kpc', shade = True)

ax.semilogy(r, densmw, label = 'mw')
# ax.semilogy(r, densmw0, label = 'mw 0')
# ax.semilogy(r, densmw1, label = 'mw 1')
# ax.semilogy(r, densmw2, label = 'mw 2')
ax.semilogy(r, denssum, label = 'mw NFW sum')
# anatomy of a plot #
plt.legend()
ax.set_title('DM Halo Mass Density')
ax.set_xlabel(r'radius ($kpc$)')
ax.set_ylabel(r'mass density ($M_{sol}/{\rm pc}^{3}$)')
ax.set_xlim((0, 10**1.5))
plt.grid()
plt.show()


# # Graph mass ##
# # NOTE profile seems shallow but m200 and r200 seem correct.
# fig, ax = plt.subplots()

# ax.semilogy(r, nfwp.mass_profile, label = 'NFW_8.8kpc')
# ax.axvline(nfwp.r_200.value, color='black', lw=0.5)
# ax.axhline(nfwp.m_200.value, color='black', lw=0.5)
# plt.legend()
# ax.set_title('DM Halo Mass')
# ax.set_xlabel(r'radius ($kpc$)')
# ax.set_ylabel(r'mass ($M_{sol}$)')

# plt.grid()
# plt.show()

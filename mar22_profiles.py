import profile as pf

import numpy as np
import matplotlib.pyplot as plt

from astropy import units as un
from astropy import constants as cnst

r_s   = 8.8  * un.kpc
v_h   = 220  * un.km / un.s
m     = 6E11 * un.Msun
rho_s = 1    * un.Msun / (un.pc)**3
r     = np.logspace(-3, 1.1, 256) * un.kpc

gp    = pf.McMillanProfile(r_s = r_s, rho_s = rho_s, v_h = v_h)

density = gp.density(r)

plt.semilogy(r, density, label = 'NFW') # semilogx
plt.legend()
plt.xlabel(r'radius ($kpc$)')
plt.ylabel(r'mass density ($M_{sol}/{\rm pc}^{3}$)')

plt.show()
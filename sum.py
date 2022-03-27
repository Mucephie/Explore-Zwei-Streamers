import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import solve_ivp 
import astropy.units as un
import astropy.constants as cnst
import profile as pf

# Defining the ODEs 
G= 4.541*10**-30
k = 1 
f_DM = 1 # for simplicity
prop = pf.McMillanProfile(r_s = 15 * un.kpc, rho_s = 0.011 * un.Msun / (un.pc)**3 )
r = 15 * un.kpc
rho = prop.density(r) # cant be greater than 0.02 for around sun
print('stream DM density: ', rho)
print('local density: ', prop.density(8.8 * un.kpc))
t = np.linspace(0, 12 * 10E7, int(10E6)) # * un.Gyr cnst.G

GG = cnst.G.to(un.pc**3 / (un.Msun * un.yr**2))
# print(GG)

def f(t, y): 
    vel, width = y 
    #make sure input and output units match

    vel = vel * un.km / un.s
    width = width * un.pc
    t = t * un.yr
    
    # unless t and width can beat out truncation error VVV this wont work
    x = (2*np.pi*GG*k*rho*f_DM)**2 * t * width # look into  units
    # should vel be delta vel?
    z = 2 * np.sqrt(width * vel) # look into units 

    fx = np.array((x.value, z.value))
    return fx


# need units for width
sol = solve_ivp(f, [0, 12 * 10E7], [0.001, 0.001], t_eval = t)
velocity = sol.y[0]
width = sol.y[1]
tt = sol.t

fig, ax = plt.subplots(1, 2)
ax1, ax2 = ax[0], ax[1]

ax1.semilogy((tt * un.yr).to(un.Gyr), velocity, label = 'velocity dispersion ')
ax1.set_xlabel('time (Gyr)')
ax1.set_ylabel('velocity')

ax2.semilogy((tt * un.yr).to(un.Gyr), width, label = 'width dispersion')
ax2.set_xlabel('time (Gyr)')
ax2.set_ylabel('width')

plt.legend()
plt.show()


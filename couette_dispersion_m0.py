import numpy as np
import matplotlib.pyplot as plt 
import scipy
from scipy import optimize
from scipy import special

nmax = 301

m = 0 # Use other script for m = 1
nnk = 101
kmax = 6.0
k = np.linspace(0,kmax,nnk) # vertical wavenumber
Omega = 1.0
a = 1.0

bess_zeros = scipy.special.jn_zeros(1, 4)

# omega valid for R = 1.0

nu0 = bess_zeros[0]
omega0 = m*Omega + 2.0 * Omega / np.sqrt((nu0**2 / k**2) + 1.0)

nu1 = bess_zeros[1]
omega1 =  m*Omega + 2.0 * Omega / np.sqrt(((nu1**2 / k**2) + 1.0))

nu2 = bess_zeros[2]
omega2 =  m*Omega +2.0 * Omega / np.sqrt(((nu2**2 / k**2) + 1.0))

nu3 = bess_zeros[3]
omega3 =  m*Omega +2.0 * Omega / np.sqrt(((nu3**2 / k**2) + 1.0))

# Negative omega


omega0m = m*Omega - 2.0 * Omega / np.sqrt((nu0**2 / k**2) + 1.0)

omega1m =  m*Omega - 2.0 * Omega / np.sqrt(((nu1**2 / k**2) + 1.0))

omega2m =  m*Omega -2.0 * Omega / np.sqrt(((nu2**2 / k**2) + 1.0))

omega3m =  m*Omega -2.0 * Omega / np.sqrt(((nu3**2 / k**2) + 1.0))

# Group speeds

cg0 = np.zeros(nnk) 
cg1 = np.zeros(nnk)  
cg2 = np.zeros(nnk)  

dk =  kmax / nnk

cg0[0] = (omega0[1]-omega0[0]) / (k[1] - k[0])
cg0[1:-2] = 0.5 / dk * ( omega0[2:-1] - omega0[0:-3] )
cg0[nnk-2] = (omega0[nnk-2]-omega0[nnk-3]) / (k[nnk-2] - k[nnk-3])
cg0[nnk-1] = (omega0[nnk-1]-omega0[nnk-2]) / (k[nnk-1] - k[nnk-2])


cg1[0] = (omega1[1]-omega1[0]) / (k[1] - k[0])
cg1[1:-2] = 0.5 / dk * ( omega1[2:-1] - omega1[0:-3] )
cg1[nnk-2] = (omega1[nnk-2]-omega1[nnk-3]) / (k[nnk-2] - k[nnk-3])
cg1[nnk-1] = (omega1[nnk-1]-omega1[nnk-2]) / (k[nnk-1] - k[nnk-2])

cg2[0] = (omega2[1]-omega2[0]) / (k[1] - k[0])
cg2[1:-2] = 0.5 / dk * ( omega2[2:-1] - omega2[0:-3] )
cg2[nnk-2] = (omega2[nnk-2]-omega2[nnk-3]) / (k[nnk-2] - k[nnk-3])
cg2[nnk-1] = (omega2[nnk-1]-omega2[nnk-2]) / (k[nnk-1] - k[nnk-2])

# Wave frequencies

plt.figure(figsize=(8,8), facecolor = 'w', edgecolor = 'k')

ax = plt.subplot(3, 1, 1)
plt.grid(True)

plt.plot(k, omega0, 'r', linewidth = 1.5, label='$\omega_1$')
plt.plot(k, omega1, 'k', linewidth = 1.5, label='$\omega_2$')
plt.plot(k, omega2, 'b',linewidth = 1.5, label='$\omega_3$')

plt.plot(k, omega0m, 'r', linestyle='--', linewidth = 1.5)
plt.plot(k, omega1m, 'k', linestyle='--',linewidth = 1.5)
plt.plot(k, omega2m, 'b', linestyle='--', linewidth = 1.5)

plt.ylim(-2,2)
plt.legend()
plt.ylabel('$\omega$') 

# Phase speed
    
ax = plt.subplot(3, 1, 2)

plt.grid(True)
plt.plot(k, omega0/k, 'r', linewidth = 1.5, label='$c_1$')
plt.plot(k, omega1/k, 'k', linewidth = 1.5, label='$c_2$')
plt.plot(k, omega2/k, 'b',linewidth = 1.5, label='$c_3$')
plt.ylim(0,0.6)
plt.ylabel('c')
plt.legend()

ax = plt.subplot(3, 1, 3)

plt.grid(True)
plt.plot(k, cg0, 'r', linewidth = 1.5, label='$c_g^1$')
plt.plot(k, cg1, 'k', linewidth = 1.5, label='$c_g^2$')
plt.plot(k, cg2, 'b',linewidth = 1.5, label='$c_g^3$')

plt.ylim(0,0.6)
plt.xlabel('$k$')
plt.ylabel('$c_g$')
plt.legend()

if (m == 0):
  plt.savefig("./vessel_dispersion_m0.eps", format = 'eps')

plt.show()


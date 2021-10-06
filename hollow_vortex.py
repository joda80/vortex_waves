import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
from scipy import special

#==============================================================================
# Kelvin's hollow vortex solution for helical waves and Kelvin's approximations;
# This solution is not discussed in the paper, but simply use Eq. (41) 
# with u_2/g_2, and note that p_1 = 0 and dP_1/dr = 0.  For u_2 and p_2 use Eqs. (122) 
# and (123), respectively, and solve for omega (which is contained in 
# g_2 = omega - m*Omega); the result is the "Full solution" below.
#
# There are no radial modes; for each m there is just one cograde and
# one retrograde mode.  The only countergrade mode (omega < 0 for all k) 
# exists for m = 1.
#==============================================================================

# Settings

R     = 1.0
Omega = 1.0
kmax = 10.0

# ====================================================
# Start script
# ====================================================

k = np.linspace(0,kmax, 4000) 

plt.figure(figsize=(8,10), facecolor = 'w', edgecolor = 'k')

for m in [1, 2, 3]:

# Full solution

  Km  = scipy.special.kv (m, k*R)
  Kpr = scipy.special.kvp(m, k*R)

  omega_co    = Omega * (m + np.sqrt(-k*R * Kpr/Km))
  omega_retro = Omega * (m - np.sqrt(-k*R * Kpr/Km)) 

# Approximations for m = 1

  N = 1.0 + k**2 * R**2 * (np.log(1.0/(k*R)) + 0.1159) # Approximation for K'/K

  omega_retro_approx = Omega * (m - np.sqrt(N))
  omega_co_approx    = Omega * (m + np.sqrt(N))

# Approximations oscillation frequency, Kelvin's Eqs. (43) and (44)

  omega_approx_retro2 = -0.5 * Omega * k**2 * R**2 * (np.log(1.0/(k*R)) + 0.1159)
  omega_approx_co2    = 2.0*Omega + 0.5 * Omega * k**2 * R**2 * (np.log(1.0/(k*R)) + 0.1159)

# =============================================================================
# Plot dispersion relations for m = 1, 2, 3 (m = 1 is solid, m > 1 is dashed)
# Legend: Black: retrograde; red: cograde;
#         Solid: m = 1; dashed: m = 2; dash-dotted: m = 3;
#         For m = 1, the asymptotic solutions are plotted in blue for small k.
# =============================================================================

  # Oscillation frequencies

  plt.subplot(3, 1, 1)
  plt.grid(True)
 
  if (m == 1):
    plt.plot(k, omega_co, color = 'r',  linewidth = 1.5)
    plt.plot(k, omega_retro, color='k', linewidth = 1.5)

  elif (m == 2):
    plt.plot(k, omega_co, color = 'r',  linewidth = 1.5, linestyle='--')
    plt.plot(k, omega_retro, color='k', linewidth = 1.5, linestyle='--')

  elif (m == 3):
    plt.plot(k, omega_co, color = 'r',  linewidth = 1.5, linestyle='-.')
    plt.plot(k, omega_retro, color='k', linewidth = 1.5, linestyle='-.')


  if (m == 1): # Approximations only for m = 1

    plt.plot(k[0:400], omega_co_approx[0:400], color='b', linewidth = 1.0)
    plt.plot(k[0:400], omega_retro_approx[0:400], color='b', linewidth = 1.0)

    plt.plot(k[0:400], omega_approx_retro2[0:400], color='b', linewidth = 1.0, linestyle='--')
    plt.plot(k[0:400], omega_approx_co2[0:400], color='b', linewidth = 1.0, linestyle='--')

  plt.xlim(0,kmax)
  plt.ylim(-4,8)

  plt.ylabel('$\omega$ (s$^{-1}$)')

# Azimuthal angular phase speed

  plt.subplot(3, 1, 2)
  plt.grid(True)

  if (m == 1):
    plt.plot(k, omega_co/m, color = 'r',  linewidth = 1.5)
    plt.plot(k, omega_retro/m, color='k', linewidth = 1.5)

  elif (m == 2):
    plt.plot(k, omega_co/m, color = 'r',  linewidth = 1.5, linestyle='--')
    plt.plot(k, omega_retro/m, color='k', linewidth = 1.5, linestyle='--')

  elif (m == 3):
    plt.plot(k, omega_co/m, color = 'r',  linewidth = 1.5, linestyle='-.')
    plt.plot(k, omega_retro/m, color='k', linewidth = 1.5, linestyle='-.')


  plt.xlim(0,kmax)
  plt.ylim(-3,5)

  plt.ylabel('$\omega/m$')

# Axial phase speed

  plt.subplot(3, 1, 3)
  plt.grid(True)

  if (m == 1):
    plt.plot(k, omega_co/k, color = 'r',  linewidth = 1.5)
    plt.plot(k, omega_retro/k, color='k', linewidth = 1.5)

  elif (m == 2):
    plt.plot(k, omega_co/k, color = 'r',  linewidth = 1.5, linestyle='--')
    plt.plot(k, omega_retro/k, color='k', linewidth = 1.5, linestyle='--')

  elif (m == 3):
    plt.plot(k, omega_co/k, color = 'r',  linewidth = 1.5, linestyle='-.')
    plt.plot(k, omega_retro/k, color='k', linewidth = 1.5, linestyle='-.')


  plt.xlim(0,kmax)
  plt.ylim(-1,4)

  plt.ylabel('$c_z$ (m s$^{-1}$)')
  plt.xlabel('k (m$^{-1}$)')

plt.show()

# The End.

import numpy as np
import matplotlib.pyplot as plt 
import scipy
import cmath
import math
from scipy import optimize
from scipy import special
from scipy import misc

# Solving Kelvin's (1880) equation 50 (nondimensional version)
# Here we define beta, find the solution of the dispersion
# relation, and then calculate omega.
#
# This script plots the dispersion relation for
# spiral modes in a Rankine vortex. The 
# root finder requires a high resolution of the
# function (nmax must be large) to resolve the peaks.
# The advantage over the built-in root finder is that
# no prior information is needed about where the root
# is located.
#
# If the roots cannot be identified, the wave frequency
# will be zero; in this case adjust nmax and nnk to
# ensure that all peaks are sufficiently resolved.
#
# Questions/comments: johannes.dahl@ttu.edu

#------------------------------------------------
# Set parameters
#------------------------------------------------

m = 1  # m >0

# Best not to change these:

kmax     =  6.0   # maximum axial wavenumber
nmax     = 100001 # Resolution of the curves whose intersections
                  # are to be found
beta_max = 12.0   # maximum eigenvalue
nnk      = 36     # For how many k

#-------------------------------------------------
# Start
#-------------------------------------------------

omega0_arr = np.zeros(nnk)
omega1_arr = np.zeros(nnk)
omega2_arr = np.zeros(nnk)
omega3_arr = np.zeros(nnk)

k_arr  = np.zeros(nnk)

c1 = 0 # Counter

plt.figure(figsize=(8,10), facecolor = 'w', edgecolor = 'k')

for k in np.linspace(0.05,kmax,nnk):
  beta0     = np.zeros(10) # Roots of dispersion relation 
  val       = np.zeros(10)   
  omega_n   = np.zeros(10)
  beta_root = np.zeros(10)

  # Initialize values of interest for complex omega

  beta_arr = np.zeros(nmax)
  lhs1d = np.zeros(shape=(nmax))
  rhs1d = np.zeros(shape=(nmax))

  counter = 0

  for mm in np.arange(nmax):

# Incrementing real argument beta  

      beta = beta_max/(np.real(nmax-1)) * np.real(mm)

      if (beta == 0):
        beta = 1.E-12

      beta_arr[mm] = beta

# Calculating omega and g

      omega = m - 2.0 * k / np.sqrt(beta**2 + k**2)

      g    = 1.0E-12 + omega  - m 

# Bessel function (1st kind) and its derivative

      bess     = scipy.special.jv (m, beta)
      bess_der = scipy.special.jvp(m, beta)

# Modified Bessel function (2nd kind) and its derivative

      mbess     = scipy.special.kv (m, k)
      mbess_der = scipy.special.kvp(m, k)

# Dispersion relation

      lhs = 1.0/beta *  bess_der/(bess+1.0E-12) 
      rhs = -1.0/ k * mbess_der/(mbess+1.0E-12) + 2.0*m / (g * beta**2) 

      lhs1d[mm] = lhs
      rhs1d[mm] = rhs

# End loop through beta

  # Finding the intersections in decreasing branches.
 
  counter = 0

  disp = lhs1d - rhs1d
  for i in np.arange(nmax-1):
    if (disp[i] > 0 and disp[i+1] <= 0 and beta_arr[i]):
      val  [counter] = lhs1d[i] # value of lhs where lhs = rhs
      beta0[counter] = beta_arr[i]
      counter = counter + 1

#---------------------------------------------------
# Plot graphical version of the solution for some k
#---------------------------------------------------

  ax = plt.subplot(4,1,1)

  if (k >=1.0 and k < 1.2 and beta0[0] > 0):
#  if (k == 5.15):

    y = np.linspace(0, beta_max, nmax)

    plt.grid(True)
    plt.plot(y, rhs1d, 'b', linewidth = 1.0) 
    plt.plot(beta0[0], val[0], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='c')
    plt.plot(beta0[1], val[1], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='r')
    plt.plot(beta0[2], val[2], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='g')
    plt.plot(beta0[3], val[3], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='b')
    plt.plot(y, lhs1d, 'k', linewidth = 1.0) 
    plt.ylim(-1,1)
    plt.ylabel('lhs (black) and rhs (blue)')
#----------------------------------------
# Calculate wave frequency 
#----------------------------------------

  omega_n = m - 2.0 * k / np.sqrt(beta0**2 + k**2)

  print(c1, k, beta0[0:4]) 

  omega0_arr[c1] = omega_n[0]
  omega1_arr[c1] = omega_n[1]
  omega2_arr[c1] = omega_n[2]
  omega3_arr[c1] = omega_n[3]

  k_arr     [c1]  = k

  c1 = c1 + 1

print(omega0_arr)
print(k_arr)

# Obtaining the group speeds (centered difference; at end points: one-sided differences)

cg0 = 0.0*omega0_arr
dk  =  kmax / nnk

cg0    [0] = (omega0_arr[1]-omega0_arr[0]) / (k_arr[1] - k_arr[0])
cg0 [1:-2] = 0.5 / dk * ( omega0_arr[2:-1] - omega0_arr[0:-3] ) 
cg0[nnk-2] = (omega0_arr[nnk-2]-omega0_arr[nnk-3]) / (k_arr[nnk-2] - k_arr[nnk-3])
cg0[nnk-1] = (omega0_arr[nnk-1]-omega0_arr[nnk-2]) / (k_arr[nnk-1] - k_arr[nnk-2])

# --

cg1 = 0.0*omega1_arr

cg1    [0] = (omega1_arr[1]-omega1_arr[0]) / (k_arr[1] - k_arr[0])
cg1 [1:-2] = 0.5 / dk * ( omega1_arr[2:-1] - omega1_arr[0:-3] )
cg1[nnk-2] = (omega1_arr[nnk-2]-omega1_arr[nnk-3]) / (k_arr[nnk-2] - k_arr[nnk-3])
cg1[nnk-1] = (omega1_arr[nnk-1]-omega1_arr[nnk-2]) / (k_arr[nnk-1] - k_arr[nnk-2])

#--

cg2 = 0.0*omega2_arr

cg2    [0] = (omega2_arr[1]-omega2_arr[0]) / (k_arr[1] - k_arr[0])
cg2 [1:-2] = 0.5 / dk * ( omega2_arr[2:-1] - omega2_arr[0:-3] )
cg2[nnk-2] = (omega2_arr[nnk-2]-omega2_arr[nnk-3]) / (k_arr[nnk-2] - k_arr[nnk-3])
cg2[nnk-1] = (omega2_arr[nnk-1]-omega2_arr[nnk-2]) / (k_arr[nnk-1] - k_arr[nnk-2])

cg3 = 0.0*omega3_arr

cg3    [0] = (omega3_arr[1]-omega3_arr[0]) / (k_arr[1] - k_arr[0])
cg3 [1:-2] = 0.5 / dk * ( omega3_arr[2:-1] - omega3_arr[0:-3] )
cg3[nnk-2] = (omega3_arr[nnk-2]-omega3_arr[nnk-3]) / (k_arr[nnk-2] - k_arr[nnk-3])
cg3[nnk-1] = (omega3_arr[nnk-1]-omega3_arr[nnk-2]) / (k_arr[nnk-1] - k_arr[nnk-2])

#--------------------------------------
# Plot wave frequency and phase speeds
#--------------------------------------

ax = plt.subplot(4,1,2)

if (m == 0):
  plt.grid(True)
  plt.plot(k_arr, omega0_arr, 'r', linewidth = 2.0, label='$\omega_0$')
  plt.plot(k_arr, omega1_arr, 'g', linewidth = 2.0, label='$\omega_1$', alpha = 0.8)
  plt.plot(k_arr, omega2_arr, 'b', linewidth = 2.0, label='$\omega_2$', alpha = 0.6)
  plt.legend()
  plt.ylabel('$\omega$')

  ax = plt.subplot(4,1,3)

  plt.grid(True)
  plt.plot(k_arr, omega0_arr/k_arr, 'r', linewidth = 2.0, label='$c_0$')
  plt.plot(k_arr, omega1_arr/k_arr, 'g', linewidth = 2.0, label='$c_1$') 
  plt.plot(k_arr, omega2_arr/k_arr, 'b', linewidth = 2.0, label='$c_2$')
  plt.xlim(0.0,6.0)

elif (m == 1 or m == -1):
  plt.grid(True)
  plt.plot(k_arr, omega0_arr, 'c', linewidth = 1.4, label='$\omega_0$')
  plt.plot(k_arr, omega1_arr, 'r', linewidth = 1.4, label='$\omega_1$')
  plt.plot(k_arr, omega2_arr, 'g', linewidth = 1.4, label='$\omega_2$')
  plt.plot(k_arr, omega3_arr, 'b', linewidth = 1.4, label='$\omega_3$')
  plt.ylim(-1.5,1.5) 
  plt.legend(loc='upper right', prop={'size': 7})

  plt.ylabel('$\omega$')
  plt.xlim(0.0,6.0)

# ---

  ax = plt.subplot(4,1,3)

  plt.grid(True)
  plt.plot(k_arr, omega0_arr/k_arr, 'c', linewidth = 1.4, label='$c_0$')
  plt.plot(k_arr, omega1_arr/k_arr, 'r', linewidth = 1.4, label='$c_1$')
  plt.plot(k_arr, omega2_arr/k_arr, 'g', linewidth = 1.4, label='$c_2$')
  plt.plot(k_arr, omega3_arr/k_arr, 'b', linewidth = 1.4, label='$c_3$')
  plt.xlim(0.0,6.0)

  plt.ylim(-0.5,2.0)
  plt.legend(loc='upper right', prop={'size': 7})
  plt.ylabel('$c$')

# --

  ax = plt.subplot(4,1,4)
  plt.grid(True)

  plt.plot(k_arr, cg0, 'c', linewidth = 1.4, label='$c_0$')
  plt.plot(k_arr, cg1, 'r', linewidth = 1.4, label='$c_1$') 
  plt.plot(k_arr, cg2, 'g', linewidth = 1.4, label='$c_2$')
  plt.plot(k_arr, cg3, 'b', linewidth = 1.4, label='$c_3$')

  plt.ylim(-0.6,0.0)
  plt.xlim(0.0,6.0)

  plt.ylabel('$c_g$')
  plt.legend(loc='lower right', prop={'size': 7})
  plt.xlabel('$k$')

if (m == 1):
  plt.savefig("./rankine_pure_m1_retro.eps", format = 'eps')

plt.show()
quit()


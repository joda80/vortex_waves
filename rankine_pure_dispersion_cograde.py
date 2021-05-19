import numpy as np
import matplotlib.pyplot as plt 
import scipy
import cmath
import math
from scipy import optimize
from scipy import special
from scipy import misc
import shapely
from shapely import geometry
from shapely.geometry import LineString, Point

# Solving Kelvin's equation 50 (nondimensional version)
# Here we define beta, find the solution of the dispersion
# relation, and then calculate omega.
#
#------------------------------------------------
# Set parameters
#------------------------------------------------

m = 2 # Limited to positive wavenumbers;

# Best not to change these:

kmax =  6.0 #  11.0
nmax = 100001
beta_max = 14.0#11.0

nnk = 29

#-------------------------------------------------
# Start
#-------------------------------------------------

Omega = 1.0

#beta = 3.8317
#r0 = 1.20051152
#u0 = -0.25655284
#jv = scipy.special.jv(1, beta*r0)
#w0 = 1.0/np.pi * u0 * 3.8317/jv
#print(w0)
#print(beta/np.pi)
#
#quit()

omega0_arr = np.zeros(nnk)
omega1_arr = np.zeros(nnk)
omega2_arr = np.zeros(nnk)

k_arr  = np.zeros(nnk)

c1 = 0 # Counter

plt.figure(figsize=(8,10), facecolor = 'w', edgecolor = 'k')

for k in np.linspace(0.0,kmax,nnk):
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

      omega = m + 2.0 * k / np.sqrt(beta**2 + k**2)

      g    = 1.0E-12 + omega  - m 

# Bessel function (1st kind) and its derivative

      bess     = scipy.special.jv (m, beta)
      bess_der = scipy.special.jvp(m, beta)

# Modified Bessel function (2nd kind) and its derivative

      mbess     = scipy.special.kv (m, k)
      mbess_der = scipy.special.kvp(m, k)

# Dispersion relation

# Little test, lhs going to zero for large k:
#      for i in np.arange(4):
#        dk = 40.0
#        mbess     = scipy.special.kv (m, dk*i+1.E-12)
#        mbess_der = scipy.special.kvp(m, dk*i+1.E-12)
#        rhs = -1.0/ (dk * i+1.E-12) * mbess_der/(mbess+1.0E-12)
#        print(i, dk*i, rhs)
#      quit()


      lhs = 1.0/beta *  bess_der/(bess+1.0E-12)
      rhs = -1.0/ k * mbess_der/(mbess+1.0E-12) + 2.0*m / (g * beta**2) 

      lhs1d[mm] = lhs
      rhs1d[mm] = rhs

# End loop through beta

  # Finding the intersections in decreasing branches.
  # First reduce the rhs so it intersects the imperfectly resolved
  # peaks.  This increses the solution just a little, but close enough
  # for our purposes.
 
  counter = 0

  for i in np.arange(nmax-1):
    if (lhs1d[i] >= rhs1d[i] and lhs1d[i+1] < rhs1d[i]):
      val  [counter] = lhs1d[i] # value of lhs where lhs = rhs
      beta0[counter] = beta_arr[i]
      counter = counter + 1

#---------------------------------------------------
# Plot graphical version of the solution for some k
#---------------------------------------------------

  ax = plt.subplot(3,1,1)

#  if ((k > 0.1 and c1 % 30 == 0)):
  if (k >= 1.0 and k < 1.1):
    y = np.linspace(0, beta_max, nmax)

    plt.grid(True)
    plt.plot(y, rhs1d, 'b', linewidth = 1.0) 
    plt.plot(beta0[0], val[0], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='r')
    plt.plot(beta0[1], val[1], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='g')
    plt.plot(beta0[2], val[2], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='b')
    print('yy', k)
    plt.plot(y, lhs1d, 'k', linewidth = 1.0) 
    plt.ylim(-1,5)
#    plt.xlabel('$\beta$')
    plt.ylabel('rhs (black) and lhs (blue)')

#----------------------------------------
# Calculate wave frequency 
#----------------------------------------

  omega_n = m +  2.0 * k / np.sqrt(beta0**2 + k**2)

  if (k == 0):
    omega_n[0] = m*Omega
    omega_n[1] = m*Omega 
    omega_n[2] = m*Omega

  print('H', c1, k, beta0, omega_n[0]) #, omega_n[0],  np.sqrt(beta0[0]**2 + k**2))


  omega0_arr[c1] = omega_n[0]
  omega1_arr[c1] = omega_n[1]
  omega2_arr[c1] = omega_n[2]

  k_arr     [c1]  = k

  c1 = c1 + 1

print(k_arr)
print()
print(omega0_arr)
print()
print(omega0_arr/k_arr)

# Obtaining the group speeds (centered difference; at end points: one-sided differences)

cg0 = 0.0*omega0_arr
dk =  kmax / nnk

cg0[0] = (omega0_arr[1]-omega0_arr[0]) / (k_arr[1] - k_arr[0])
cg0[1:-2] = 0.5 / dk * ( omega0_arr[2:-1] - omega0_arr[0:-3] ) 
cg0[nnk-2] = (omega0_arr[nnk-2]-omega0_arr[nnk-3]) / (k_arr[nnk-2] - k_arr[nnk-3])
cg0[nnk-1] = (omega0_arr[nnk-1]-omega0_arr[nnk-2]) / (k_arr[nnk-1] - k_arr[nnk-2])

# --

cg1 = 0.0*omega1_arr

cg1[0] = (omega1_arr[1]-omega1_arr[0]) / (k_arr[1] - k_arr[0])
cg1[1:-2] = 0.5 / dk * ( omega1_arr[2:-1] - omega1_arr[0:-3] )
cg1[nnk-2] = (omega1_arr[nnk-2]-omega1_arr[nnk-3]) / (k_arr[nnk-2] - k_arr[nnk-3])
cg1[nnk-1] = (omega1_arr[nnk-1]-omega1_arr[nnk-2]) / (k_arr[nnk-1] - k_arr[nnk-2])

#--

cg2 = 0.0*omega2_arr

cg2[0] = (omega2_arr[1]-omega2_arr[0]) / (k_arr[1] - k_arr[0])
cg2[1:-2] = 0.5 / dk * ( omega2_arr[2:-1] - omega2_arr[0:-3] )
cg2[nnk-2] = (omega2_arr[nnk-2]-omega2_arr[nnk-3]) / (k_arr[nnk-2] - k_arr[nnk-3])
cg2[nnk-1] = (omega2_arr[nnk-1]-omega2_arr[nnk-2]) / (k_arr[nnk-1] - k_arr[nnk-2])

#--------------------------------------
# Plot wave frequency and phase speeds
#--------------------------------------

ax = plt.subplot(3,1,2)

if (m == 0):
  plt.grid(True)
  plt.plot(k_arr, omega0_arr, 'r', linewidth = 2.0, label='$\omega_1$')
  plt.plot(k_arr, omega1_arr, 'g', linewidth = 2.0, label='$\omega_2$')
  plt.plot(k_arr, omega2_arr, 'b', linewidth = 2.0, label='$\omega_3$')
  plt.legend()
#plt.xlabel('$k$')
  plt.xlim(0,6)
  plt.ylabel('$\omega$')

  ax = plt.subplot(3,1,3)

  plt.grid(True)
  plt.plot(k_arr, omega0_arr/k_arr, 'r', linewidth = 2.0, label='$c_1$')
  plt.plot(k_arr, omega1_arr/k_arr, 'g', linewidth = 2.0, label='$c_2$')
  plt.plot(k_arr, omega2_arr/k_arr, 'b', linewidth = 2.0, label='$c_3$')
  ax = plt.subplot(3,1,3)


elif (m == 1 or m == -1):
  plt.grid(True)
  plt.plot(k_arr, omega0_arr, 'r', linewidth = 2.0, label='$\omega_1$')
  plt.plot(k_arr, omega1_arr, 'g', linewidth = 2.0, label='$\omega_2$', alpha = 0.8)
  plt.plot(k_arr, omega2_arr, 'b', linewidth = 2.0, label='$\omega_3$', alpha = 0.6)
  plt.legend()
#plt.xlabel('$k$')
  plt.ylabel('$\omega$')

  ax = plt.subplot(3,1,3)

  plt.grid(True)
  plt.plot(k_arr, omega0_arr/k_arr, 'r', linewidth = 2.0, label='$c_1$')
  plt.plot(k_arr, omega1_arr/k_arr, 'g', linewidth = 2.0, label='$c_2$', alpha = 0.8)
  plt.plot(k_arr, omega2_arr/k_arr, 'b', linewidth = 2.0, label='$c_3$', alpha = 0.6)
  
  plt.plot(k_arr, cg0, 'r', linewidth = 2.0, linestyle='--')
  plt.plot(k_arr, cg1, 'g', linewidth = 2.0, linestyle='--')
  plt.plot(k_arr, cg2, 'b', linewidth = 2.0, linestyle='--')
  plt.ylim(0,1.5)


plt.legend()
plt.xlabel('$k$')
plt.ylabel('$c$')

if (m == 0):
  plt.savefig("./rankine_pure_dispersion_m0.eps", format = 'eps')
if (m == 1):
  plt.savefig("./rankine_pure_dispersion_m1.eps", format = 'eps')

plt.show()
quit()



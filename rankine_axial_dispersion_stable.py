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
# 
# Corrected version of the solver for stable Kelvin modes
# in a swirling Rankinw vortex (Loiseleux scenario).
#
#------------------------------------------------
# Set parameters
#------------------------------------------------

m = 1 # Limited to positive wavenumbers

# Best not to change these:

kmax = 6.0
nmax = 90001
beta_max = 13.0 # Needs to cover all three roots
S = 1.0
nnk = 31 #11# 31 # For how many k

kmax = 2.2
nmax = 40000
nnk = 2

W = 1.0
Omega = 1.0

#-------------------------------------------------
# Start
#-------------------------------------------------

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

omega0m_arr = np.zeros(nnk)
omega1m_arr = np.zeros(nnk)
omega2m_arr = np.zeros(nnk)

k_arr  = np.zeros(nnk)

c1 = 0 # Counter

plt.figure(figsize=(8,10), facecolor = 'w', edgecolor = 'k')

for k in np.linspace(0.01,kmax,nnk):
  beta0     = np.zeros(10) # Roots of dispersion relation 
  val       = np.zeros(10)   
  omega_n   = np.zeros(10)
  beta_root = np.zeros(10)

  beta0m     = np.zeros(10) # Roots of dispersion relation 
  valm       = np.zeros(10)
  beta_rootm = np.zeros(10)


  print('Beginning of the loop')

  # Initialize values of interest for complex omega

  beta_arr = np.zeros(nmax)
  lhs1d = np.zeros(shape=(nmax))
  rhs1d = np.zeros(shape=(nmax))

  lhs1dm = np.zeros(shape=(nmax))
  rhs1dm = np.zeros(shape=(nmax))

  counter = 0
  counterm = 0

  for mm in np.arange(nmax):

# Incrementing real argument beta  

      beta = beta_max/(np.real(nmax-1)) * np.real(mm)

      if (beta == 0):
        beta = 1.E-12

      beta_arr[mm] = beta

# Calculating omega and g
     
      ktr = k  # Doppler-shifted frequency because of W

      omega = m*S + ktr +  2.0 * k / np.sqrt(beta**2 + k**2)
 
      omegam = m*S + ktr -  2.0 * k / np.sqrt(beta**2 + k**2)
      
      g    = 1.0E-12 + omega -ktr  - m*S 
      gm = 1.0E-12 + omegam -ktr  - m*S

# Bessel function (1st kind) and its derivative

      bess     = scipy.special.jv (m, beta)
      bess_der = scipy.special.jvp(m, beta)

# Modified Bessel function (2nd kind) and its derivative

      mbess     = scipy.special.kv (m, k)
      mbess_der = scipy.special.kvp(m, k)

# Dispersion relation

      lhs = (g + k)**2 * beta *  bess_der/(bess+1.0E-12)
      rhs = -g**2 * beta**2 / k * mbess_der/(mbess+1.0E-12) + 2.0* m * S * (g + k)**2 / g

      lhsm = (gm + k)**2 * beta *  bess_der/(bess+1.0E-12)
      rhsm = -gm**2 * beta**2 / k * mbess_der/(mbess+1.0E-12) + 2.0* m * S * (gm + k)**2 / gm

      lhs1d[mm] = lhs
      rhs1d[mm] = rhs

      lhs1dm[mm] = lhsm
      rhs1dm[mm] = rhsm

# End loop through beta

  # Finding the intersections in decreasing branches.
  # First reduce the rhs so it intersects the imperfectly resolved
  # peaks.  This increses the solution just a little, but close enough
  # for our purposes.
 
  counter = 0

  for i in np.arange(nmax-1):
    if (lhs1d[i] >= rhs1d[i] and lhs1d[i+1] <= rhs1d[i+1] and beta_arr[i] > 1.E-4):
      val  [counter] = rhs1d[i] # value of lhs where lhs = rhs 
      beta0[counter] = beta_arr[i]
      counter = counter + 1

  counterm = 0

  for i in np.arange(nmax-1):
    if (lhs1dm[i] >= rhs1dm[i] and lhs1dm[i+1] <= rhs1dm[i+1] and beta_arr[i] > 1.E-4):
      valm  [counterm] = rhs1dm[i] # value of lhs where lhs = rhs 
      beta0m[counterm] = beta_arr[i]
      counterm = counterm + 1

#---------------------------------------------------
# Plot graphical version of the solution for some k
#---------------------------------------------------

  ax = plt.subplot(2,1,1)

#  if ((k >= 3.73 and k<=3.74) or k == 3.68): #> 0.1 and c1 % 30 == 0)):
#  if ((k >= 3.73 and k<=3.74 or k == 1.008) or k == 2.505): #> 0.1 and c1 % 30 == 0)):
  if(k >= 1.0 and k<=1.2): 
    y = np.linspace(0, beta_max, nmax)

    plt.grid(True)
    plt.plot(y, rhs1d, 'b', linewidth = 1.0)
    if (m == 0):
      plt.plot(y, rhs1dm, 'b', linestyle='--', linewidth = 2.0)
    elif (m == 1):
      plt.plot(y, rhs1dm, 'b', linestyle='--', linewidth = 1.0)

    plt.plot(beta0[0], val[0], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='r')
    plt.plot(beta0[1], val[1], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='r')
    plt.plot(beta0[2], val[2], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='r')

    plt.plot(beta0m[0], valm[0], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='g')
    plt.plot(beta0m[1], valm[1], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='g')
    plt.plot(beta0m[2], valm[2], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='g')

    plt.plot(y, lhs1d, 'k', linewidth = 1.0) 
    plt.plot(y, lhs1dm, 'k', linestyle='--', linewidth = 1.0)

    if (m == 0):
      plt.ylim(-20, 30)
    elif (m == 1):
      plt.ylim(-10,30)
#    plt.xlabel('$\beta$')
    plt.ylabel('rhs (blue) and lhs (black)')
  #  plt.show()
  #  quit()
#----------------------------------------
# Calculate wave frequency 
#----------------------------------------

  # Intrinsic frequency makes the most sense

  omega_n = m*S + k + 2.0 * k / np.sqrt(beta0**2 + k**2) - m*S  -k

  omega_nm = m*S + k - 2.0 * k / np.sqrt(beta0m**2 + k**2) -m*S -k

#  print('H', c1, k, beta0[0], omega_n[0], omega_nm[0], beta0m[0]) #, omega_n[0],  np.sqrt(beta0[0]**2 + k**2))
  print('k: ', k)
  print('beta and INTRINSIC frequency')
  print(beta0[0], omega_n[0])
  print(beta0[1], omega_n[1])
  print(beta0[2], omega_n[2])
  print('')
  print(beta0m[0], omega_nm[0])
  print(beta0m[1], omega_nm[1])
  print(beta0m[2], omega_nm[2])
  print('--------------------------')


  omega0_arr[c1] = omega_n[0]
  omega1_arr[c1] = omega_n[1]
  omega2_arr[c1] = omega_n[2]

  omega0m_arr[c1] = omega_nm[0]
  omega1m_arr[c1] = omega_nm[1]
  omega2m_arr[c1] = omega_nm[2]

 # elif (m == 1):
 #   omega0_arr[c1] = m + 2.0
 #   omega1_arr[c1] = omega_n[0]
 #   omega2_arr[c1] = omega_n[1]

  k_arr     [c1]  = k

  c1 = c1 + 1

print(omega0_arr)
print(k_arr)

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

#--
cg0m = 0.0*omega0_arr

cg0m[0] = (omega0m_arr[1]-omega0m_arr[0]) / (k_arr[1] - k_arr[0])
cg0m[1:-2] = 0.5 / dk * ( omega0m_arr[2:-1] - omega0m_arr[0:-3] )
cg0m[nnk-2] = (omega0m_arr[nnk-2]-omega0m_arr[nnk-3]) / (k_arr[nnk-2] - k_arr[nnk-3])
cg0m[nnk-1] = (omega0m_arr[nnk-1]-omega0m_arr[nnk-2]) / (k_arr[nnk-1] - k_arr[nnk-2])

# --

cg1m = 0.0*omega1_arr

cg1m[0] = (omega1m_arr[1]-omega1m_arr[0]) / (k_arr[1] - k_arr[0])
cg1m[1:-2] = 0.5 / dk * ( omega1m_arr[2:-1] - omega1m_arr[0:-3] )
cg1m[nnk-2] = (omega1m_arr[nnk-2]-omega1m_arr[nnk-3]) / (k_arr[nnk-2] - k_arr[nnk-3])
cg1m[nnk-1] = (omega1m_arr[nnk-1]-omega1m_arr[nnk-2]) / (k_arr[nnk-1] - k_arr[nnk-2])

#--

cg2m = 0.0*omega2_arr

cg2m[0] = (omega2m_arr[1]-omega2m_arr[0]) / (k_arr[1] - k_arr[0])
cg2m[1:-2] = 0.5 / dk * ( omega2m_arr[2:-1] - omega2m_arr[0:-3] )
cg2m[nnk-2] = (omega2m_arr[nnk-2]-omega2m_arr[nnk-3]) / (k_arr[nnk-2] - k_arr[nnk-3])
cg2m[nnk-1] = (omega2m_arr[nnk-1]-omega2m_arr[nnk-2]) / (k_arr[nnk-1] - k_arr[nnk-2])

# Read unstable phase speed

infile2 = './dispersion_mp1_S1_longwave.txt'

[k1, omega_r1, omega_i1, beta_unst, omega_intr]  = np.transpose(np.loadtxt(infile2, skiprows=1))

# Use intrinsic frequencies

omega_r1 = omega_r1 - m*S -k1 

# Asymptotic longwave solution, intrinsic

omega_r1_as = m*S  + k1 - k1**2/(2.0*S) - m*S -k1

#--------------------------------------
# Plot wave frequency and phase speeds
#--------------------------------------

ax = plt.subplot(3,1,2)

if (m == 0):
  plt.grid(True)
  plt.plot(k_arr, omega0_arr, 'r', linewidth = 1.5, label='$g_1$')
  plt.plot(k_arr, omega1_arr, 'g', linewidth = 1.5, label='$g_2$')
  plt.plot(k_arr, omega2_arr, 'b', linewidth = 1.5, label='$g_3$')

  plt.plot(k_arr, omega0m_arr, 'r', linestyle='--', linewidth = 1.5)
  plt.plot(k_arr, omega1m_arr, 'g', linestyle='--', linewidth = 1.5)
  plt.plot(k_arr, omega2m_arr, 'b', linestyle='--', linewidth = 1.5)
  plt.ylim(-2.0, 2.0)

  plt.legend()
#plt.xlabel('$k$')
  plt.ylabel('$\omega$')

  ax = plt.subplot(3,1,3)

  plt.grid(True)
  plt.plot(k_arr, omega0_arr/k_arr, 'r', linewidth = 1.5, label='$c_1$')
  plt.plot(k_arr, omega1_arr/k_arr, 'g', linewidth = 1.5, label='$c_2$')
  plt.plot(k_arr, omega2_arr/k_arr, 'b', linewidth = 1.5, label='$c_3$')


  plt.plot(k_arr, omega0m_arr/k_arr, 'r', linestyle='--', linewidth = 1.5)
  plt.plot(k_arr, omega1m_arr/k_arr, 'g', linestyle='--', linewidth = 1.5)
  plt.plot(k_arr, omega2m_arr/k_arr, 'b', linestyle='--', linewidth = 1.5)

  plt.ylim(-1.0, 1.0)

#  plt.plot(k_arr, cg0, 'r', linewidth = 2.0, linestyle='--')
#  plt.plot(k_arr, cg1, 'g', linewidth = 2.0, linestyle='--')
#  plt.plot(k_arr, cg2, 'b', linewidth = 2.0, linestyle='--')

elif (m == 1):
  plt.grid(True)

  plt.grid(True)
  plt.plot(k_arr, omega0_arr, 'r', linewidth = 1.5, label='$g_1$')
  plt.plot(k_arr, omega1_arr, 'k', linewidth = 1.5, label='$g_2$')
#  plt.plot(k_arr, omega2_arr, 'b', linewidth = 1.5, label='$g_3$')

  plt.plot(k_arr, omega0m_arr, 'k', linestyle='-.', linewidth = 1.5,  label='$g_0$')
  plt.plot(k_arr, omega1m_arr, 'r', linestyle='--', linewidth = 1.5, label='$g_1$')
  plt.plot(k_arr, omega2m_arr, 'k', linestyle='--', linewidth = 1.5, label='$g_2$')

  # KH mode

  plt.plot(k1, omega_r1, 'b', linestyle='-', linewidth = 2.5, label='$g_{unst}$')

#  plt.plot(k1, omega_r1_as, 'k', linestyle='--', linewidth = 2.0 )

  plt.legend()
  plt.xlim(0,7)
  plt.ylim(-2.2,2.0)
  plt.xlabel('$k (m^{-1})$')
  plt.ylabel('$g (s^{-1})$')
 
#  ax = plt.subplot(3,1,3)
#  plt.grid(True)
#  plt.plot(k_arr, omega0_arr/k_arr, 'r', linewidth = 1.5, label='$c_1$')
#  plt.plot(k_arr, omega1_arr/k_arr, 'g', linewidth = 1.5, label='$c_2$')
#  plt.plot(k_arr, omega2_arr/k_arr, 'b', linewidth = 1.5, label='$c_3$')

#  plt.plot(k_arr, omega0m_arr/k_arr, 'r', linestyle='--', linewidth = 1.5)
#  plt.plot(k_arr, omega1m_arr/k_arr, 'g', linestyle='--', linewidth = 1.5)
#  plt.plot(k_arr, omega2m_arr/k_arr, 'b', linestyle='--', linewidth = 1.5)

#  plt.ylim (-2,2.5)

#plt.legend()
#plt.xlabel('$k (s^{-1})$')
#plt.ylabel('$c$')

if (m == 0):
  plt.savefig("./kelvin_axial_dispersion_m0.eps", format = 'eps')
elif (m == 1):
  plt.savefig("./kelvin_axial_dispersion_m1.eps", format = 'eps')

plt.show()
quit()



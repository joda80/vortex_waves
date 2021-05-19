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

#------------------------------------------------
# Set parameters
#------------------------------------------------

m = 1 # For m = 0 use other script
R = 1.0
Omega = 1.0
kmax =  6.0 #  11.0
nmax = 70001
beta_max = 15.0
nnk =  26 # For how many k

#-------------------------------------------------
# Start
#-------------------------------------------------

omega0_arr = np.zeros(nnk)
omega1_arr = np.zeros(nnk)
omega2_arr = np.zeros(nnk)

omega0m_arr = np.zeros(nnk)
omega1m_arr = np.zeros(nnk)
omega2m_arr = np.zeros(nnk)

c1 = 0 # Counter

plt.figure(figsize=(7,10), facecolor = 'w', edgecolor = 'k')

counter = 0
counterm = 0

lhs = np.zeros(nmax)
rhs = np.zeros(nmax)
rhsm = np.zeros(nmax)


beta = np.linspace(0, beta_max, nmax)

k_arr = np.linspace(0,kmax, nnk)

for k in k_arr:

#  if k == 0:
#    k = 1.E-12
 
  g = 2.0*Omega / (np.sqrt( beta**2/k**2 + 1.0  ))


  bess     = scipy.special.jv (m,R*beta)
  bess_der = scipy.special.jvp(m, R*beta)

  lhs = bess_der/(bess+1.0E-12)
  rhs = 2.0*m*Omega/(g*R*beta + 1.0E-12)

  omega = m - 2.0 * k / np.sqrt(beta**2 + k**2)
  g    = 1.0E-12 + omega  - m
  rhsm = 2.0*m*Omega/(g*R*beta + 1.0E-12)

  counter = 0

  val = np.zeros(10)
  valm = np.zeros(10)

  beta0 = np.zeros(10)
  beta0m = np.zeros(10)


  for i in np.arange(nmax-1):
    if (lhs[i] > rhs[i] and lhs[i+1] <= rhs[i]):
      val  [counter] = lhs[i] # value of lhs where lhs = rhs
      beta0[counter] = beta[i]
      counter = counter + 1

# neg. omega
  counterm = 0
  for i in np.arange(nmax-1):
    if (lhs[i] > rhsm[i] and lhs[i+1] <= rhsm[i]):
      valm  [counterm] = lhs[i] # value of lhs where lhs = rhs
      beta0m[counterm] = beta[i]
      counterm = counterm + 1

  print('k   beta   beta_retro', k, beta0, beta0m)

#---------------------------------------------------
# Plot graphical version of the solution for some k
#---------------------------------------------------

  ax = plt.subplot(4,1,1)

  if (k == 1.2):
    rhs1d = np.zeros(nmax) + rhs
    rhsm1d = np.zeros(nmax) + rhsm

    plt.grid(True)
    plt.plot(beta, rhs1d, 'b', linewidth = 1.0) 
    plt.plot(beta, rhsm1d, 'b', linestyle='--', linewidth = 1.0)

    plt.plot(beta0[0], val[0], 'o',  markersize=5, markerfacecolor='w',
           markeredgewidth=1.0, markeredgecolor='r')
    plt.plot(beta0[1], val[1], 'o',  markersize=5, markerfacecolor='w',
           markeredgewidth=1.0, markeredgecolor='g')
    plt.plot(beta0[2], val[2], 'o',  markersize=5, markerfacecolor='w',
               markeredgewidth=1.0, markeredgecolor='b')

    plt.plot(beta0m[0], valm[0], 'o',  markersize=5, markerfacecolor='w',
               markeredgewidth=1.0, markeredgecolor='r', linestyle='dashed')
    plt.plot(beta0m[1], valm[1], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='g')
    plt.plot(beta0m[2], valm[2], 'o',  markersize=5, markerfacecolor='w',
             markeredgewidth=1.0, markeredgecolor='b')
    plt.plot(beta, lhs, 'k', linewidth = 1.0) 
    plt.ylim(-2,8)
    plt.xlim(0,beta_max)
#    plt.xlabel('$\beta$')
    plt.ylabel('rhs (blue) and lhs (black)')

#----------------------------------------
# Calculate wave frequency 
#----------------------------------------

  omega_n = m *Omega +  2.0 * k / np.sqrt(beta0**2 + k**2)
  omega_nm = m *Omega -  2.0 * k / np.sqrt(beta0m**2 + k**2) 

  if (k == 0):
    omega_n[0] = m *Omega
    omega_n[1] = m *Omega
    omega_n[2] = m *Omega

    omega_nm[0] = m *Omega
    omega_nm[1] = m *Omega
    omega_nm[2] = m *Omega

  if (m == 1):
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

# Retrograde modes

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

#--------------------------------------
# Plot wave frequency and phase speeds
#--------------------------------------

# Frequency

ax = plt.subplot(4,1,2)

plt.grid(True)
plt.plot(k_arr, omega0_arr, 'r', linewidth = 1.5, label='$\omega_1$')
plt.plot(k_arr, omega1_arr, 'g', linewidth = 1.5, label='$\omega_2$')
plt.plot(k_arr, omega2_arr, 'b', linewidth = 1.5, label='$\omega_3$')

plt.plot(k_arr, omega0m_arr, 'r', linestyle = '--', linewidth = 1.5)
plt.plot(k_arr, omega1m_arr, 'g', linestyle = '--', linewidth = 1.5)
plt.plot(k_arr, omega2m_arr, 'b', linestyle = '--', linewidth = 1.5)
plt.xlim(0,kmax)
plt.legend()
#plt.xlabel('$k$')
plt.ylabel('$\omega$')

# Phase speed

ax = plt.subplot(4,1,3)

plt.grid(True)
plt.plot(k_arr, omega0_arr/k_arr, 'r', linewidth = 1.5, label='$c_1$')
plt.plot(k_arr, omega1_arr/k_arr, 'g', linewidth = 1.5, label='$c_2$')
plt.plot(k_arr, omega2_arr/k_arr, 'b', linewidth = 1.5, label='$c_3$')

plt.plot(k_arr, omega0m_arr/k_arr, 'r', linestyle='--', linewidth = 1.5)
plt.plot(k_arr, omega1m_arr/k_arr, 'g', linestyle='--', linewidth = 1.5)
plt.plot(k_arr, omega2m_arr/k_arr, 'b', linestyle='--', linewidth = 1.5)
plt.ylabel('$c$')


plt.legend(loc='upper right')
plt.ylim(0, kmax)
plt.ylim(-0.5,2.0)
plt.xlim(0,kmax)
# Group speed

ax = plt.subplot(4,1,4)

plt.grid(True)

plt.plot(k_arr, cg0, 'r', linewidth = 1.5, label='$c_1$')
plt.plot(k_arr, cg1, 'g', linewidth = 1.5)
plt.plot(k_arr, cg2, 'b', linewidth = 1.5)

plt.plot(k_arr, cg0m, 'r', linewidth = 1.5, linestyle='--', label='$c^g_1$')
plt.plot(k_arr, cg1m, 'g', linewidth = 1.5, linestyle='--')
plt.plot(k_arr, cg2m, 'b', linewidth = 1.5, linestyle='--')
plt.xlim(0,kmax)

plt.legend()
plt.xlabel('$k$')
plt.ylabel('$c_g$')
plt.xlim(0,kmax)

if (m == 1):
  plt.savefig("./vessel_dispersion_m1.eps", format = 'eps')

plt.show()
quit()



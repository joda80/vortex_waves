import numpy as np
import matplotlib.pyplot as plt 
import scipy
from scipy import optimize
from scipy import special

# Plot of the axisymmetric Kelvin waves in the
# meridional (r,z) plane. 

#--------------------------------------------
# Settings
#--------------------------------------------

icase = 1 # Bounded solid-bocy rotation (Couette flow)
#icase = 2 # Rankine vortex
nmax = 301

#--------------------------------------------

u0 = 1.0 # Were we prescribe u0
t  = 0.0

#-------------------------------------------------------------------------
# Case 1: Couette modes
#-------------------------------------------------------------------------

if icase == 1:
  rmax = 1.0
  zmax = 10.0
  rmat = 0.6 # Material boundary
  bd_scale = 0.2 
  k = 1.5
  beta = 3.83170597/rmax
  omega = 0.73 
  Omega = 1.0
  w0 = beta/k
  rmaxp = rmax + 1.0

  r = np.linspace(-rmaxp,rmaxp,nmax)  # Technically, r*beta
  z = np.linspace(0,zmax,nmax)

  R,Z = np.meshgrid(r,z)

  u =  np.zeros(shape=(nmax,nmax))
  v = np.zeros(shape=(nmax,nmax))
  w = np.zeros(shape=(nmax,nmax))
  J1 = scipy.special.jv (1, beta*R)

  g = omega

# Couette, fundamental mode

  ilim = np.zeros(2)
  ilim[0] = rmax
  ilim[1] = -rmax

  zd = np.linspace(0,zmax, nmax)
  rd = rmat - bd_scale * w0/g * k / beta * scipy.special.jvp (0, beta*rmat) * np.cos(k*zd-omega*t)

  for i in np.arange(nmax):

    if r[i] >= 0 and r[i]<= ilim[0]:
      u[:, i] =   u0 * scipy.special.jv (1, beta*R[:, i]) * np.sin(k*Z[:, i]- omega*t)
      v[:,i]  = - w0 * 2.0*Omega*k/(beta*omega) * scipy.special.jv (1, beta*R[:, i])*np.cos(k*Z[:, i]- omega*t)    
      w[:, i] =   w0 * scipy.special.jv (0, beta*R[:, i]) * np.cos(k*Z[:, i]- omega*t)

    if r[i] >= ilim[1] and r[i] <= 0:
      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.sin(k*Z[:, i]- omega*t)
      v[:,i]  = -w0 * 2.0*Omega*k/(beta*omega) * scipy.special.jv (1, beta*R[:, i])*np.cos(k*Z[:, i]- omega*t)
      w[:, i] =  w0 * scipy.special.jv (0, beta*R[:, i]) * np.cos(k*Z[:, i]- omega*t)    

    if r[i] >= 0 and r[i] > ilim[0]:
      u[:, i] = 0.0 
      v[:, i] = 0.0
      w[:, i] = 0.0
    if r[i] <=0 and r[i] < ilim[1]:
      u[:, i] = 0.0
      v[:, i] = 0.0
      w[:, i] = 0.0


#---------------------------------------------------------
# Plot Couette fundamental mode
#---------------------------------------------------------

  plt.figure(figsize=(10,10), facecolor = 'w', edgecolor = 'k')

  # Left panel

  plt.subplot(2,2,1)

  istart = 0
  iend   = nmax
  istride = 8

  zstart = 0
  zend = nmax
  
  plt.vlines(rmax, 0, zmax, linewidths=2, colors='k')
  plt.vlines(-rmax, 0, zmax, linewidths=2, colors='k')

  plt.quiver(R[zstart:zend:istride, istart:iend:istride], Z[zstart:zend:istride, istart:iend:istride],
  u[zstart:zend:istride, istart:iend:istride], w[zstart:zend:istride, istart:iend:istride], color = '0.3',scale = 40.0)

  plt.plot(rd,zd, color='r', linewidth=2)
  plt.plot(-rd,zd, color='r', linewidth=2)
  plt.xlim(-rmaxp, rmaxp)
  plt.ylim(0,zmax)
  plt.ylabel('z (m)')


  # Right panel

  plt.subplot (2, 2,2)

  plt.plot(rd,zd, color='r', linewidth=2)
  plt.plot(-rd,zd, color='r', linewidth=2)

  plt.vlines(rmax, 0, zmax, linewidths=2, colors='k')
  plt.vlines(-rmax, 0, zmax, linewidths=2, colors='k')

  plt.streamplot(r[istart:iend], z,
  u[zstart:zend, istart:iend], w[zstart:zend, istart:iend],
  density=2.6,arrowsize=0.5, linewidth=0.3, color='k')
  plt.xlim(-rmaxp, rmaxp)
  plt.ylim(0,zmax)

# ---------------------------
# Couette, third radial mode
# ---------------------------

  rmax = 10.0
  zmax = 10.0
  r1 = 3.0 # Material boundary
  r2 = 6.0
  r3 = 9.0
  bd_scale = 0.5 
  k = 1.5
  beta = 10.17346814/rmax
  omega = 0.3 # dispersion relation
  w0 = beta/k
  rmaxp = rmax + 2.0

  r = np.linspace(-rmaxp,rmaxp,nmax)  # Technically, r*beta
  z = np.linspace(0,zmax,nmax)

  R,Z = np.meshgrid(r,z)

  u =  np.zeros(shape=(nmax,nmax)) #  
  w = np.zeros(shape=(nmax,nmax))
  J1 = scipy.special.jv (1, beta*R)

  g = omega

# Couette, fundamental mode

  ilim = np.zeros(2)
  ilim[0] = rmax
  ilim[1] = -rmax

  zd = np.linspace(0,zmax, nmax)
  rd1 = r1 - bd_scale * w0/g * k / beta * scipy.special.jvp (0, beta*r1) * np.cos(k*Z-omega*t)
  rd2 = r2 - bd_scale * w0/g * k / beta * scipy.special.jvp (0, beta*r2) * np.cos(k*Z-omega*t)
  rd3 = r3 - bd_scale * w0/g * k / beta * scipy.special.jvp (0, beta*r3) * np.cos(k*Z-omega*t)


  for i in np.arange(nmax):

    if r[i] >= 0 and r[i]<= ilim[0]:
      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.sin(k*Z[:, i]- omega*t)
      w[:, i] = w0 * scipy.special.jv (0, beta*R[:, i]) * np.cos(k*Z[:, i]- omega*t)

    if r[i] >= ilim[1] and r[i] <= 0:
      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.sin(k*Z[:, i]- omega*t)
      w[:, i] =  w0 * scipy.special.jv (0, beta*R[:, i]) * np.cos(k*Z[:, i]- omega*t)

    if r[i] >= 0 and r[i] > ilim[0]:
      u[:, i] = 0.0
      w[:, i] = 0.0
    if r[i] <=0 and r[i] < ilim[1]:
      u[:, i] = 0.0
      w[:, i] =  0.0

#--------------------------------------
# Plot Couette third mode
#--------------------------------------

  istride = 6

  # Left panel

  plt.subplot (2, 2,3)

  plt.vlines(rmax, 0, zmax, linewidths=2, colors='k')
  plt.vlines(-rmax, 0, zmax, linewidths=2, colors='k')

  plt.quiver(R[zstart:zend:istride, istart:iend:istride], Z[zstart:zend:istride, istart:iend:istride],
  u[zstart:zend:istride, istart:iend:istride], w[zstart:zend:istride, istart:iend:istride], color = '0.3',scale = 20.0)

  plt.plot(rd1,zd, color='r', linewidth=2)
  plt.plot(-rd1,zd, color='r', linewidth=2)

  plt.plot(rd2,zd, color='r', linewidth=2)
  plt.plot(-rd2,zd, color='r', linewidth=2)

  plt.plot(rd3,zd, color='r', linewidth=2)
  plt.plot(-rd3,zd, color='r', linewidth=2)

  plt.xlim(-rmaxp, rmaxp)
  plt.ylim(0,zmax)
  plt.ylabel('z (m)')
  plt.xlabel('r (m)')

  # Right panel

  plt.subplot (2, 2,4)

  plt.plot(rd1,zd, color='r', linewidth=2)
  plt.plot(-rd1,zd, color='r', linewidth=2)

  plt.plot(rd2,zd, color='r', linewidth=2)
  plt.plot(-rd2,zd, color='r', linewidth=2)

  plt.plot(rd3,zd, color='r', linewidth=2)
  plt.plot(-rd3,zd, color='r', linewidth=2)

  plt.vlines(rmax, 0, zmax, linewidths=2, colors='k')
  plt.vlines(-rmax, 0, zmax, linewidths=2, colors='k')

  plt.streamplot(r[istart:iend], z,
  u[zstart:zend, istart:iend], w[zstart:zend, istart:iend],
  density=2.6,arrowsize=0.5, linewidth=0.3, color='k')
  plt.xlim(-rmaxp, rmaxp)
  plt.ylim(0,zmax)
  plt.xlabel('r (m)')

  plt.savefig("./vessel_structure_m0.eps", format = 'eps')

  plt.show()
  quit()

#---------------------------------------------------------------------
# Case 2: Rankine vortex (fundamental and third modes)
#---------------------------------------------------------------------

elif icase == 2:
  
  RMW = 1.0
  rmax = RMW + 1.0

# Allowable combinations

#  k = 5:    beta = 3.25625, 
#  k = 2:    beta = 2.87925
#  k = 1:    beta = 2.6509, 
#  k = 0.5:  beta = 2.513
#  k = 0.05: beta = 2.408  
  
  k = 0.05  # Pick k here

  if (k == 0.05):
    beta = 2.408
    zmax = 180.0
    ssc = 1200
    istride = 10
    bd_scale = 0.075


  elif (k == 1.0):
    beta = 2.6509
    zmax = 10.0
    ssc = 40
    bd_scale = 0.3
    istride = 6

  elif (k == 5.0):
    beta =  3.25625 # From kelvin_dispersion.py
    zmax = 2.0
    ssc = 20
    bd_scale = 1.0 
    istride = 6

  B = 1.0

  beta = beta/RMW
  u0 = 1.0
  w0 = beta/k

  r = np.linspace(-rmax,rmax,nmax)
  z = np.linspace(0,zmax,nmax)

  R,Z = np.meshgrid(r,z)

  u = np.zeros(shape=(nmax,nmax))
  w = np.zeros(shape=(nmax,nmax))
  uout = np.zeros(shape=(nmax,nmax))

  ilim = np.zeros(2)
  ilim[0] = RMW 
  ilim[1] = -RMW 
 
  alpha = beta * RMW
  Omega = 1
  g = 2.0 * Omega / (np.sqrt(alpha**2 / (RMW*k)**2 + 1.0))

  omega = g

  Jm_R  = scipy.special.jv (0, beta*RMW)
  Km_R  = scipy.special.kv (0, beta*RMW)
  Kpr_R = scipy.special.kvp(0, beta*RMW)

  zd = np.linspace(0,zmax, nmax)
  rd = RMW + bd_scale * w0/g * Kpr_R/Km_R * Jm_R * np.cos(k*zd-omega*t)

  for i in np.arange(nmax):

#    if ((i % 50 == 0)):
#      print(i)

    if r[i] >= 0 and r[i]< ilim[0]:
      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.sin(k*Z[:, i]- omega*t)
      w[:, i] =  w0 * scipy.special.jv (0, beta*R[:, i]) * np.cos(k*Z[:, i]- omega*t)

    if r[i] > ilim[1] and r[i] <= 0:
      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.sin(k*Z[:, i]- omega*t)
      w[:, i] = w0 * scipy.special.jv (0, beta*R[:, i]) * np.cos(k*Z[:, i]- omega*t)    

  J0R = scipy.special.jv (0, beta*RMW)
  K0R = scipy.special.kv (0, k*RMW)
  B   = w0 * J0R/K0R
  
  for i in np.arange(nmax):

#    if ((i % 50 == 0)):
#      print(i)

    if r[i] >= 0 and r[i] >= ilim[0]: 
      u[:, i] = B *  scipy.special.kvp (0, k*R[:, i]) * np.sin(k*Z[:, i]- omega*t)
      w[:, i] =  B * scipy.special.kv (0, k*R[:, i]) * np.cos(k*Z[:, i]- omega*t)
    if r[i] <=0 and r[i] <= ilim[1]: # r < 0, so we need to flip signs
      u[:, i] =  -B  * scipy.special.kvp (0, -k*R[:, i]) * np.sin(k*Z[:, i]- omega*t)
      w[:, i] =  B * scipy.special.kv (0, -k*R[:, i]) * np.cos(k*Z[:, i]- omega*t)

  J0R = scipy.special.jv (0, beta*RMW)
  K0R = scipy.special.kv (0, k*RMW)

#------------------------------------------------------
# Plot Rankine radial modes (each k separately)
#------------------------------------------------------

  istart = 0
  iend   = nmax-1 

  zstart = 0
  zend = nmax

  plt.figure(figsize=(10, 5), facecolor = 'w', edgecolor = 'k')
  
  plt.subplot (1, 2,1)
  plt.vlines(ilim, 0, zmax, linestyle='--',linewidths=0.5, colors='k')

  if (k == 0.05):

    plt.quiver(R[zstart:zend:istride, istart:iend:istride], Z[zstart:zend:istride, istart:iend:istride],
    u[zstart:zend:istride, istart:iend:istride], w[zstart:zend:istride, istart:iend:istride], 
    width = 0.005, headwidth = 3.5, color = '0.3',scale = ssc)

    plt.plot(rd,zd, color='r', linewidth=2)
    plt.plot(-rd,zd, color='r', linewidth=2)
   
    plt.xlabel('r (m)')
    plt.ylabel('z (m)')
    plt.xlim(-2, 2)

    plt.subplot (1, 2,2)

    plt.plot(rd,zd, color='r', linewidth=2)
    plt.plot(-rd,zd, color='r', linewidth=2)

    plt.streamplot(r[istart:iend], z,
    u[zstart:zend, istart:iend], w[zstart:zend, istart:iend],
    density=1.3,arrowsize=0.6, linewidth=0.3, color='k')

  elif (k > 0.05):
 
    plt.quiver(R[zstart:zend:istride, istart:iend:istride], Z[zstart:zend:istride, istart:iend:istride],
    u[zstart:zend:istride, istart:iend:istride], w[zstart:zend:istride, istart:iend:istride], 
    color = '0.3',scale = ssc)
  
    plt.plot(rd,zd, color='r', linewidth=2)
    plt.plot(-rd,zd, color='r', linewidth=2)

    if (k == 1.0):
      plt.xlabel('r (m)')
    plt.ylabel('z (m)')
    plt.xlim(-2, 2)

    plt.subplot (1, 2,2)
    plt.streamplot(r[istart:iend], z,   
    u[zstart:zend, istart:iend], w[zstart:zend, istart:iend],
    density=1.5,arrowsize=0.6, linewidth=0.3, color='k')

    plt.plot(rd,zd, color='r', linewidth=2)
    plt.plot(-rd,zd, color='r', linewidth=2)

  plt.vlines(ilim, 0,zmax, linestyle='--',linewidths=0.5, colors='k')
  plt.xlabel('r (m)')
  plt.xlim(-2, 2)
 
  if (k == 5.0): 
    plt.savefig("./meridional_structure_k500.eps", format = 'eps')

  elif (k == 1.0):
    plt.savefig("./meridional_structure_k100.eps", format = 'eps')

  elif (k == 0.05):
    plt.savefig("./meridional_structure_k005.eps", format = 'eps')

  plt.show()
  quit()

# All done.

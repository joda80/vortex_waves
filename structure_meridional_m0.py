import numpy as np
import matplotlib.pyplot as plt 
import scipy
from scipy import optimize
from scipy import special

#--------------------------------------------
# Settings
#--------------------------------------------

icase = 1
nmax = 301

#--------------------------------------------

u0 = 1.0
t = 0.0

#-------------------------------------------------------------------------
# Case 1: Couette modes
#-------------------------------------------------------------------------

if icase == 1:
  rmax = 1.0
  zmax = 10.0
  rmat = 0.6 # Material boundary
  bd_scale = 0.2 
  k2 = 1.5
  k = k2
  beta = 3.83170597/rmax
  omega = 1.1 # Correct is 0.73, but this practically only affects the amplitude of the boundary
  Omega = 1.0
  w0 = beta/k2
  d0 = 0.25
  rmaxp = rmax + 1.0

  r = np.linspace(-rmaxp,rmaxp,nmax)  # Technically, r*beta
  z = np.linspace(0,zmax,nmax)

  R,Z = np.meshgrid(r,z)

  u =  np.zeros(shape=(nmax,nmax)) #  u0 * scipy.special.jv (1, R) * np.sin(k2*Z)
  v = np.zeros(shape=(nmax,nmax))
  w = np.zeros(shape=(nmax,nmax))
  J1 = scipy.special.jv (1, beta*R)

  g = omega

# Couette, fundamental mode
  ilim = np.zeros(2)
  ilim[0] = rmax
  ilim[1] = -rmax

  zd = np.linspace(0,zmax, nmax)
  rd = rmat - bd_scale * w0/g * k / beta * scipy.special.jvp (0, beta*rmat) * np.cos(k2*zd-omega*t)

  for i in np.arange(nmax):

    if r[i] >= 0 and r[i]<= ilim[0]:
      u[:, i] =   u0 * scipy.special.jv (1, beta*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      v[:,i]  = - w0 * 2.0*Omega*k2/(beta*omega) * scipy.special.jv (1, beta*R[:, i])*np.cos(k2*Z[:, i]- omega*t)    
      w[:, i] =   w0 * scipy.special.jv (0, beta*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)

    if r[i] >= ilim[1] and r[i] <= 0:
      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      v[:,i]  = -w0 * 2.0*Omega*k2/(beta*omega) * scipy.special.jv (1, beta*R[:, i])*np.cos(k2*Z[:, i]- omega*t)
      w[:, i] =  w0 * scipy.special.jv (0, beta*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)    

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

#  plt.contour(R,Z,v)

  plt.plot(rd,zd, color='r', linewidth=2)
  plt.plot(-rd,zd, color='r', linewidth=2)
  plt.xlim(-rmaxp, rmaxp)
  plt.ylim(0,zmax)
  plt.ylabel('z (m)')


#  plt.contour(1.0/beta * R[::istride, istart:iend:istride], Z[::istride, istart:iend:istride],
#  d[::istride, istart:iend:istride],
#  linewidths = 2, linestyles='solid', colors = 'r',  levels = [ilim[1]/beta+0.5, ilim[0]/beta-0.5])

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
  k2 = 1.5
  k = k2
  beta = 10.17346814/rmax
  omega = 0.3 # dispersion relation
  w0 = beta/k2
  d0 = 0.25
  rmaxp = rmax + 2.0

  r = np.linspace(-rmaxp,rmaxp,nmax)  # Technically, r*beta
  z = np.linspace(0,zmax,nmax)

  R,Z = np.meshgrid(r,z)

  u =  np.zeros(shape=(nmax,nmax)) #  u0 * scipy.special.jv (1, R) * np.sin(k2*Z)
  w = np.zeros(shape=(nmax,nmax))
  J1 = scipy.special.jv (1, beta*R)

  g = omega

# Couette, fundamental mode
  ilim = np.zeros(2)
  ilim[0] = rmax
  ilim[1] = -rmax

  zd = np.linspace(0,zmax, nmax)
  rd1 = r1 - bd_scale * w0/g * k / beta * scipy.special.jvp (0, beta*r1) * np.cos(k2*Z-omega*t)
  rd2 = r2 - bd_scale * w0/g * k / beta * scipy.special.jvp (0, beta*r2) * np.cos(k2*Z-omega*t)
  rd3 = r3 - bd_scale * w0/g * k / beta * scipy.special.jvp (0, beta*r3) * np.cos(k2*Z-omega*t)


  for i in np.arange(nmax):

    if r[i] >= 0 and r[i]<= ilim[0]:
      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      w[:, i] = w0 * scipy.special.jv (0, beta*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)

    if r[i] >= ilim[1] and r[i] <= 0:
      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      w[:, i] =  w0 * scipy.special.jv (0, beta*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)

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

#  plt.contour(1.0/beta * R[::istride, istart:iend:istride], Z[::istride, istart:iend:istride],
#  d[::istride, istart:iend:istride],
#  linewidths = 2, linestyles='solid', colors = 'r',  levels = [ilim[1]/beta+0.5, ilim[0]/beta-0.5])

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

#  bess0_zeros = scipy.special.jn_zeros(0, 3)

#  b0 = bess0_zeros[0]
  
  RMW = 1.0
  rmax = RMW + 1.0

# Allowable combinations


#  k = 5:    beta = 3.25625, B = -11.5
#  k = 2:    beta = 2.87925
#  k = 1:    beta = 2.6509, B = -0.86
#  k = 0.5:  beta = 2.513
#  k = 0.05: beta = 2.408  B = -0.52
#  k2 = 1.0
#k2 = 0.05
  
  k2 = 0.05

  if (k2 == 0.05):
    beta = 2.408
    zmax = 180.0
    ssc = 1200
    istride = 10
    bd_scale = 0.075


  elif (k2 == 1.0):
    beta = 2.6509
    zmax = 10.0
    ssc = 40
    bd_scale = 0.3
    istride = 6

  elif (k2 == 5.0):
    beta =  3.25625 # From kelvin_dispersion.py
    zmax = 2.0
    ssc = 20
    bd_scale = 1.0 
    istride = 6

  B = 1.0

  beta = beta/RMW
  u0 = 1.0
  w0 = beta/k2

  r = np.linspace(-rmax,rmax,nmax)
  z = np.linspace(0,zmax,nmax)

  R,Z = np.meshgrid(r,z)

  u = np.zeros(shape=(nmax,nmax)) #u0 * scipy.special.jv (1, beta*R) * np.sin(k2*Z)
  w = np.zeros(shape=(nmax,nmax))
  uout = np.zeros(shape=(nmax,nmax))

 # bess1_zeros = scipy.special.jn_zeros(1, 3)
  ilim = np.zeros(2)
  ilim[0] = RMW 
  ilim[1] = -RMW 
 

# Proper formulation of the boundary displacement (equivalent to the ine implemented)
 
  alpha = beta * RMW
  Omega = 1
  g = 2.0 * Omega / (np.sqrt(alpha**2 / (RMW*k2)**2 + 1.0))

  omega = g

  Jm_R  = scipy.special.jv (0, beta*RMW)
  Km_R  = scipy.special.kv (0, beta*RMW)
  Kpr_R = scipy.special.kvp(0, beta*RMW)

  zd = np.linspace(0,zmax, nmax)
  rd = RMW + bd_scale * w0/g * Kpr_R/Km_R * Jm_R * np.cos(k2*zd-omega*t)

#  plt.plot(z,rd)
#  plt.contour(R,Z,R)
#  plt.plot(rd, zd)
#  plt.show()
#  quit()

  print('Loop 1')
  for i in np.arange(nmax):

    if ((i % 50 == 0)):
      print(i)

    if r[i] >= 0 and r[i]< ilim[0]:
      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      w[:, i] =  w0 * scipy.special.jv (0, beta*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)

    if r[i] > ilim[1] and r[i] <= 0:
      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      w[:, i] = w0 * scipy.special.jv (0, beta*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)    

  print('Loop 2')

  J0R = scipy.special.jv (0, beta*RMW)
  K0R = scipy.special.kv (0, k2*RMW)
  print(B)
  B = w0 * J0R/K0R
  print(B)
  print()
  
  for i in np.arange(nmax):

    if ((i % 50 == 0)):
      print(i)

    if r[i] >= 0 and r[i] >= ilim[0]: 
      u[:, i] = B *  scipy.special.kvp (0, k2*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      w[:, i] =  B * scipy.special.kv (0, k2*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)
    if r[i] <=0 and r[i] <= ilim[1]: # r < 0, so we need to flip signs
      u[:, i] =  -B  * scipy.special.kvp (0, -k2*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      w[:, i] =  B * scipy.special.kv (0, -k2*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)

  J0R = scipy.special.jv (0, beta*RMW)
  K0R = scipy.special.kv (0, k2*RMW)

#------------------------------------------------------
# Plot Rankine radial modes (each k separately)
#------------------------------------------------------

  istart = 0#  5
  iend   = nmax-1 # nmax-5

  zstart = 0
  zend = nmax

  plt.figure(figsize=(10, 5), facecolor = 'w', edgecolor = 'k')
  
  plt.subplot (1, 2,1)
  plt.vlines(ilim, 0, zmax, linestyle='--',linewidths=0.5, colors='k')

  if (k2 == 0.05):

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

  elif (k2 > 0.05):
 
    plt.quiver(R[zstart:zend:istride, istart:iend:istride], Z[zstart:zend:istride, istart:iend:istride],
    u[zstart:zend:istride, istart:iend:istride], w[zstart:zend:istride, istart:iend:istride], 
    color = '0.3',scale = ssc)
  
    plt.plot(rd,zd, color='r', linewidth=2)
    plt.plot(-rd,zd, color='r', linewidth=2)

    if (k2 == 1.0):
      plt.xlabel('r (m)')
    plt.ylabel('z (m)')
    plt.xlim(-2, 2)

    plt.subplot (1, 2,2)
    plt.streamplot(r[istart:iend], z,   
    u[zstart:zend, istart:iend], w[zstart:zend, istart:iend],
    density=1.5,arrowsize=0.6, linewidth=0.3, color='k')

    plt.plot(rd,zd, color='r', linewidth=2)
    plt.plot(-rd,zd, color='r', linewidth=2)

#    plt.contourf(R[:, istart:iend], Z[:, istart:iend],
#    w[zstart:zend, istart:iend], levels=np.linspace(-3,3,21))
#    plt.colorbar()

  plt.vlines(ilim, 0,zmax, linestyle='--',linewidths=0.5, colors='k')
  plt.xlabel('r (m)')
  plt.xlim(-2, 2)
 
  if (k2 == 5.0): 
    plt.savefig("./meridional_structure_k500.eps", format = 'eps')

  elif (k2 == 1.0):
    plt.savefig("./meridional_structure_k100.eps", format = 'eps')

  elif (k2 == 0.05):
    plt.savefig("./meridional_structure_k005.eps", format = 'eps')

  plt.show()
  quit()

#-----------------------------------------------------------------

elif icase == 4:

  RMW = 1.0
  k2 = 1.0
  beta = 5.642/RMW
  w0 = beta/RMW
  d0 = 1.0

  r = np.linspace(-1.5,1.5,nmax)
  z = np.linspace(0,10,nmax)

  R,Z = np.meshgrid(r,z)

  ilim = np.zeros(2)
  ilim[0] = RMW 
  ilim[1] = -RMW 


  u = np.zeros(shape=(nmax,nmax)) # u0 * scipy.special.jv (1, beta*R) * np.sin(k2*Z)
  w = np.zeros(shape=(nmax,nmax))

  for i in np.arange(nmax):

    d = R - d0/omega * scipy.special.jv (1, beta*R) * np.cos(k2*Z-omega*t)

    # Right side, core

    if r[i] >= 0 and r[i]<= ilim[0]:
      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)
      w[:, i] = -w0 * scipy.special.jv (0, beta*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)

    # Left side, core

    if r[i] >= ilim[1] and r[i] <= 0:
      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)
      w[:, i] =  -w0 * scipy.special.jv (0, beta*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)

    # Right side, outer region

    B = -0.52 # k = 1
    if r[i] >= 0 and r[i] > ilim[0]:
      u[:, i] =  -B * k2 * scipy.special.kvp (0, k2*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)
      w[:, i] =  -B * k2 * scipy.special.kv (0, k2*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)

    # Left side, outer region

    if r[i] <=0 and r[i] < ilim[1]:
      u[:, i] = B * k2 * scipy.special.kvp (0, -k2*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)
      w[:, i] = -B * k2 * scipy.special.kv (0, -k2*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
 
# Match solutions at the boundary visually
#   
  if (12==1):
    B = -0.52 
    for i in np.arange(nmax):
#      u[:, i] = u0 * scipy.special.jv (1, beta*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)
#      w[:, i] =  -B*k2* scipy.special.kvp (0, k2*R[:, i]) * np.cos(k2*Z[:, i]- omega*t)
      u[:, i] = -w0 * scipy.special.jv (0, beta*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      w[:, i] = B * k2 * scipy.special.kv (0, k2*R[:, i]) * np.sin(k2*Z[:, i]- omega*t)

#
    plt.plot(r[150:nmax], u[10, 150:nmax], 'r')
    plt.plot(r[150:nmax], w[10, 150:nmax], 'b')
    plt.plot(r[150:nmax], 0* w[10, 150:nmax], 'k', '--', linewidth = 0.7)
    plt.vlines(ilim[0], -2, 2)
    plt.xlabel('Radial distance')
    plt.ylabel('Inner (red) and outer (blue) solutions for w')
    plt.ylim(-2,2)

    plt.show()
    quit()

# Swirling Rankine flow

if icase == 5:

  r = np.linspace(-10,10,nmax)
  z = np.linspace(0,10,nmax)

  R,Z = np.meshgrid(r,z)

  k = 1.0
  W = 1.0
  g = (omega - k*W)

  u = u0 * scipy.special.jv (1, R) * np.sin(k2*Z)
  w = np.zeros(shape=(nmax,nmax))

  bess1_zeros = scipy.special.jn_zeros(1, 3)

  print(bess1_zeros)
  ilim = np.zeros(2)
  ilim[0] =  bess1_zeros[0]
  ilim[1] = -bess1_zeros[0]

  for i in np.arange(nmax):

    d = R - d0/omega * scipy.special.jv (1, R) * np.cos(k2*Z-omega*t)

    if r[i] >= 0 and r[i]<= ilim[0]:
      u[:, i] = u0 * scipy.special.jv (1, R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      w[:, i] =  W + w0 * scipy.special.jv (0, R[:, i]) * np.cos(k2*Z[:, i]- omega*t)

    if r[i] >= ilim[1] and r[i] <= 0:
      u[:, i] = u0 * scipy.special.jv (1, R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      w[:, i] = W + w0 * scipy.special.jv (0, R[:, i]) * np.cos(k2*Z[:, i]- omega*t)

  for i in np.arange(nmax):

    # u should be zero because that's where we put te boundary.  But since the
    # modified 2nd kind Bessel function drops off so rapidly, it is practically zero (about 0.02).

    if r[i] >= 0 and r[i] > ilim[0]:
      u[:, i] =  0.0 # -1.0 * scipy.special.kvp (0, R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      w[:, i] = (W*k/(g + 3.0E-14) + 1.0) * 1.E-12 * (-k2) * scipy.special.kv(0, R[:, i]) * np.cos(k2*Z[:, i]- omega*t)
    if r[i] <=0 and r[i] < ilim[1]:
      u[:, i] = 0.0 # -1.0 * scipy.special.kvp (0, -R[:, i]) * np.sin(k2*Z[:, i]- omega*t)
      w[:, i] = (W*k/(g + 3.0E-14) + 1.0) * 1.E-12 * (-k2) * scipy.special.kv(0, -R[:, i]) * np.cos(k2*Z[:, i] - omega*t)

# Match solutions visually, using Loieleux matching condition, u+ = (u-)*(WK/g + 1)
# k = 1; W = 1; still not working.

#    u[:, i] =  W+ w0 * scipy.special.jv (0, R[:, i]) * np.cos(k2*Z[:, i]- omega*t) 
#    w[:, i] = (W*k/(g + 3.0E-14) + 1.0) * 15.0E-13 * ( (k2) * scipy.special.kv(0, R[:, i]) * np.cos(k2*Z[:, i]- omega*t))

#  plt.plot(r[150:nmax], u[150, 150:nmax], 'r')
#  plt.plot(r[150:nmax], w[150, 150:nmax], 'b')
#  plt.plot(r[150:nmax], 0* w[150, 150:nmax], 'k', linewidth = 0.7)
#  plt.vlines(ilim[0], -2, 2)
#  plt.xlabel('Radial distance')
#  plt.ylabel('Inner (red) and outer (blue) solutions for w')
#  plt.ylim(-2,2)

#  plt.show()
#  quit()

#----------------------
#if icase == 1:
#  dr = 20.0/np.real(nmax)
#elif icase == 2:
#  dr = 24.0/np.real(nmax)

#dudr = np.zeros(shape=(nmax,nmax))
#dudr[:, 1:-2] = 0.5 / dr * (u[:,2:-1] - u[:,0:-3])

#dz = 5.0/np.real(nmax)
#dwdz = np.zeros(shape=(nmax,nmax))
#dwdz[1:-2,:] = 0.5 / dz * (w[2:-1,:] - w[0:-3,:])
#div = beta*u/(R+1.E-12) + beta*dudr + dwdz
#-------------------------

if (icase == 1):

  istart = 0
  iend   = nmax
  istride = 6

  zstart = 0
  zend = nmax
  fig, ax = plt.subplots()
  plt.vlines(ilim/beta, 0, 10, linewidths=2, colors='k')
  plt.quiver(1.0/beta*R[zstart:zend:istride, istart:iend:istride], Z[zstart:zend:istride, istart:iend:istride], 
  u[zstart:zend:istride, istart:iend:istride], w[zstart:zend:istride, istart:iend:istride], color = '0.3',scale = 30.0) 
  
  plt.contour(1.0/beta * R[::istride, istart:iend:istride], Z[::istride, istart:iend:istride],
  d[::istride, istart:iend:istride],
  linewidths = 2, linestyles='solid', colors = 'r',  levels = [ilim[1]/beta+0.5, ilim[0]/beta-0.5])


  plt.xlabel('Radial distance')
  plt.ylabel('Height')

elif (icase == 2):
  istart = 0#  5
  iend   = nmax-1 # nmax-5

  zstart = 0
  zend = nmax

  plt.figure(figsize=(5, 10), facecolor = 'w', edgecolor = 'k')
  #fig, ax = plt.subplots(2, 1)
  plt.subplot (2, 1,1)
  plt.vlines(ilim, 0, zmax, linestyle='--',linewidths=0.5, colors='k')
  plt.quiver(R[zstart:zend:istride, istart:iend:istride], Z[zstart:zend:istride, istart:iend:istride],
  u[zstart:zend:istride, istart:iend:istride], w[zstart:zend:istride, istart:iend:istride], color = '0.3',scale = ssc)

  plt.contour(R[:, istart:iend], Z[:, istart:iend],
  d[:, istart:iend],
  linewidths = 2, linestyles='solid', colors = 'r',  levels =[-RMW,RMW])
  plt.xlabel('Radial distance')
  plt.ylabel('Height')
  plt.xlim(-rmax, rmax)


  plt.subplot (2, 1,2)
  plt.streamplot(r[istart:iend], z,
  u[zstart:zend, istart:iend], w[zstart:zend, istart:iend],
  density=1.5,arrowsize=0.6, linewidth=0.3, color='k')

  plt.vlines(ilim, 0,zmax, linestyle='--',linewidths=0.5, colors='k')
  plt.xlabel('Radial distance')
  plt.ylabel('Height')
  plt.xlim(-rmax, rmax)

elif (icase == 4):

  istart = 0#  5
  iend   = nmax-1 # nmax-5

  zstart = 0
  zend = nmax
  zmax = 5.0

  istride = 10
  ssc = 60.0
 

  plt.figure(figsize=(5, 10), facecolor = 'w', edgecolor = 'k')
  plt.subplot (2, 1,1)
  plt.vlines(ilim, 0, 10, linestyle='--',linewidths=1.5, colors='k')
  plt.quiver(R[zstart:zend:istride, istart:iend:istride], Z[zstart:zend:istride, istart:iend:istride],
  u[zstart:zend:istride, istart:iend:istride], w[zstart:zend:istride, istart:iend:istride], color = '0.3',scale = ssc)

  levs = [-10,-4.2, 4.2, 10]
  plt.contour(R[zstart:zend:istride, istart:iend:istride], Z[zstart:zend:istride, istart:iend:istride],
  d[zstart:zend:istride, istart:iend:istride],
  linewidths = 2, linestyles='solid', colors = 'r',  levels = levs)
  plt.xlabel('Radial distance')
  plt.ylabel('Height')

  plt.subplot (2, 1,2)
  plt.streamplot(r[istart:iend], z,
  u[zstart:zend, istart:iend], w[zstart:zend, istart:iend],
  density=1.5,arrowsize=0.6, linewidth=0.3, color='k')

  plt.vlines(ilim, 0, 10, linestyle='--',linewidths=1.5, colors='k')
  plt.xlabel('Radial distance')
  plt.ylabel('Height')

if (icase == 3):
  plt.savefig("./pics/kelvin_structure_k%02d.eps" %(k2), format = 'eps')
elif (icase == 4):
  plt.savefig("./pics/kelvin_structure_mode3_k%02d.eps" %(k2), format = 'eps')

plt.show()



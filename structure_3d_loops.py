import netCDF4
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt 
import scipy
from scipy import optimize
from scipy import special
import plotly
import plotly.graph_objects as go
import numpy as np
from matplotlib.colors import LightSource
from matplotlib import cm
import pylab

from numpy import sin, cos, pi
from skimage import measure
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#=========================================================================
# BOUNDED VORTEX (Couette)
#=========================================================================

icase = 4

if icase == 1: # Bounded
  nmax = 101
  R = 1.0
  xmax = R + 0.5 
  ymax = R + 0.5 
  zmin = 0.0 
  zmax =  20.0

  p0 = 100000.0
  rho = 1.0
  k = 1.0
  ntmax =  100
  w0 = 1.0
  Omega = 1.0
  iretro = 0

# Bounded flow

#  m = 0
#  alpha = 3.8317
# alpha = 7.0156

  m = 1
  alpha = 4.5590875
  dr = 0.1

#  alpha = 7.7675125 
#  alpha = 10.935375 
#  dr = 0.2

# Retrograde modes

#  iretro = 1
#  m = 1
#  k = 1.2
#  alpha = 2.87207143
#  alpha = 6.09514286
#  alpha =9.26721429
#  dr = 0.2

#  m = 2
#  alpha = 5.5794375 
#  alpha  = 8.86875
#  alpha = 12.074625
#  dr = 0.1

#  m = 4
#  alpha =  7.828875
#  alpha = 11.306625
#  alpha = 14.61525
#  dr = 0.1

# -------------------
# Done with settings
# -------------------

  beta = alpha/R
  
  for nt in np.arange(0, ntmax):

    X,Y,Z = np.mgrid[(-nmax/2.0):(nmax/2.0):1, (-nmax/2.0):(nmax/2.0):1, 0:nmax:1]
    X = X * 2.0 * xmax/nmax
    Y = Y * 2.0 * ymax/nmax
    Z = Z * zmax/ nmax

    XX = X * np.cos(np.pi) + Y * np.sin(np.pi)
    YY = X * np.sin(np.pi) + Y * np.cos(np.pi)

    THETA =  np.pi + np.arctan2(YY,XX)

    RR = np.sqrt(X**2 + Y**2 + 0*Z**2)

# Calculate frequencies and velocities

    g = 2.0 * Omega / (np.sqrt(alpha**2 / (R*k)**2 + 1.0))

    if (iretro == 1):
      g = -2.0 * Omega / (np.sqrt(alpha**2 / (R*k)**2 + 1.0))
  
    omega = g + m*Omega
    print(g, omega)

    tmax = 2.0*pi/omega
    t   =  tmax/np.real(ntmax) * nt
 
    Jm  = scipy.special.jv (m, beta* RR)
    Jpr = scipy.special.jvp(m, beta*RR)

    u = -w0 * g * (g*beta*Jpr - 2.0*m*Omega/RR * Jm) / (k * (4.0 * Omega**2 - g**2)) * np.sin(k*Z + m*THETA - omega*t)
    V = Omega * RR

    vv = Omega*RR + w0 * g * (2.0*Omega * beta * Jpr - m*g/RR * Jm) / (k * (4.0 * Omega**2 - g**2)) * np.cos(k*Z + m*THETA - omega*t)
    v =  w0 * g * (2.0*Omega * beta * Jpr - m*g/RR * Jm) / (k * (4.0 * Omega**2 - g**2)) * np.cos(k*Z + m*THETA - omega*t)

    for i in np.arange(nmax):
      for j in np.arange(nmax):
        for kk in np.arange(nmax):
          if (RR[i,j,kk] > R):
            u[i,j,kk] = 0.0
            v[i,j,kk] = 0.0
            vv[i,j,kk] = 0.0
 # Pressure

    fac = 100.0
    ppert = fac*(g*rho/k)*w0*Jm * np.cos(k*Z + m*THETA - omega*t)
    pbase = p0 + rho/2.0 * Omega**2 * (RR**2 - R**2)
    p = pbase + ppert

    iso_val = np.min(np.min(p)) + 1.0
    verts, faces, dummy, dummy = measure.marching_cubes(p, iso_val, spacing=(0.1, 0.1, 0.1))

    print('Plotting...') 
    fig = plt.figure(figsize=(5,10), facecolor = 'w', edgecolor = 'k')
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=5, azim=45.0)
    ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2], antialiased=True) 

    if (nt < 10):
      tstr = '00'+str(nt)
    elif (nt >= 10 and nt < 100):
      tstr = '0'+str(nt)
    elif (nt >= 100):
      tstr = str(nt)     

    filename = "couette_pure_m1_"+tstr+".png"

    plt.savefig("./loops/"+filename) 
    print('Saved figure: ', filename)
#    plt.show()
    plt.close(fig)

  quit()

#========================================================================
# RANKINE VORTEX
#========================================================================

elif (icase == 2): # Rankine vortex

  nmax = 101
  R = 1.0
  xmax = R + 0.5
  ymax = R + 0.5
  zmin = 0.0
  zmax = 20.0
  ntmax = 100
  z = 0.0
  w0 = 1.0
  d0 = 0.2
  rho = 1.0
  p0 = 100000.0
  Omega = 1.0
  ampl = 1.0
  iunst = 0 # default--don't touch
  iretro = 0 # --"--

# Rankine

#  m = 0
#  z = np.pi # 0.0
##  k = 1.0
#  alpha = 2.650875
#  alpha = 5.642
#  alpha = 8.73325

#  m = 1
#  iretro = 1

#  k = 1.0
#  alpha = 3.957875 # fundamental
#  bd_scale = 0.8

#  alpha = 7.0915   # 1st radial mode
#  bd_scale = 0.8

# Unstable, k = 1
#  k     = 1.0
#  alpha = 9.021722615282124 # Unstable mode!
                           
#  omega = 1.7796622824320223 
#  W = 1.0
#  iunst = 1

# Unstable, k = 1.1

  k = 1.1
  m = 1
  alpha = 8.469
  omega = 1.851 
  W = 1.0
  iunst = 1
  bd_scale = 0.1

# Unstable, k = 3

#  k = 3.0
#  omega = 2.691033845651748 
#  alpha = 3.4656815907481584
#  W = 1.0
#  iunst = 1

#  k = 2.0
#  alpha = 4.11925 
#  alpha = 7.201375

# Negative intrinsic frequency 

#  iunst = 1 # = 1 for stable mode in axial ject scenario
#  iretro = 1
#  m = 1
# Countergrade

#  k = 0.05
#  alpha = 0.086  # longwave limit displacement mode

#  k = 1.07
#  alpha= 1.17312 # structureless / displacement countergrade mode 
#  bd_scale = 0.3
#  ampl = 0.5

# Retrograde modes

#  alpha=4.01688
#  bd_scale = 0.8

#  alpha=7.11444
#  alpha=10.24056

#  W = 1.0
#  k = 1.2 
#  omega = -0.287 + W + k
#  alpha = 8.329
#  bd_scale = 0.1

#  m = 2
#  k = 1.0
#  alpha = 5.2045
#  alpha = 8.4623
#  alpha =11.653425

# -------------------
# Done with settings
# -------------------

  beta = alpha/R

  for nt in np.arange(0, ntmax):

    X,Y,Z = np.mgrid[(-nmax/2.0):(nmax/2.0):1, (-nmax/2.0):(nmax/2.0):1, 0:nmax:1]
    X = X * 2.0 * xmax/nmax
    Y = Y * 2.0 * ymax/nmax
    Z = Z * zmax/ nmax

    XX = X * np.cos(np.pi) + Y * np.sin(np.pi)
    YY = X * np.sin(np.pi) + Y * np.cos(np.pi)

    THETA =  np.pi + np.arctan2(YY,XX)

    RR = np.sqrt(X**2 + Y**2)

  # Intrinsic freq

    g = 2.0 * Omega / (np.sqrt(alpha**2 / (R*k)**2 + 1.0))

  # Modify for retrograde and unstable modes
 
    if (iretro == 1):
      g = -2.0 * Omega / (np.sqrt(alpha**2 / (R*k)**2 + 1.0))

    if (iunst == 1):
      g = omega - Omega*m - W*k

  # Lab frame frequency

    omega = g + m*Omega + W*k

    Jm  = scipy.special.jv (m, beta* RR)
    Jpr = scipy.special.jvp(m, beta*RR)

    Km  = scipy.special.kv (m, k*RR)
    Kpr = scipy.special.kvp(m, k*RR)

    tmax = 2.0*pi/omega
    t   =  tmax/np.real(ntmax) * nt

    u = w0 * g * (g*beta*Jpr - 2.0*m*Omega/RR * Jm)     / (k * (4.0 * Omega**2 - g**2)) * (-np.sin(k*Z + m*THETA - omega*t))
    v = w0 * g * (2.0*Omega * beta * Jpr - m*g/RR * Jm) / (k * (4.0 * Omega**2 - g**2)) *   np.cos(k*Z + m*THETA - omega*t)
    vv = 0.8*Omega * RR + v 
 
    if (iretro == 1):
      vv = 2.5*Omega * RR + v
    if (iunst == 1):
      vv = 1.5*Omega * RR + v

    w = w0 * Jm * np.cos(k*Z + m*THETA - omega*t)

    Jpr_R = scipy.special.jvp(m, beta*R)
    Jm_R  = scipy.special.jv (m, beta* R) 
    Kpr_R = scipy.special.kvp(m, k*R)
    Km_R  = scipy.special.kv (m, k*R)

  # Matching dispacement 

    w1 = - k * w0/(g * beta**2) * (g*beta*Jpr_R - 2.0*m*Omega/R * Jm_R) / Kpr_R
    print(w1)

  # unstable
    if (iunst == 1):
   
      g1 = g
      g2 = g + W*k                          # g should be g1 but still works out 
      w1 = - k * w0 *g2/(g1**2 * beta**2) * (g1*beta*Jpr_R - 2.0*m*Omega/R * Jm_R) / Kpr_R
    #g = g + W*k # g reassigned to represent intrinsic freq outside of rising core, i.e., g2

  # w ~ p matched
   
    if (iunst == 0):
      w1 = w0 * Jm_R/Km_R
    elif (iunst == 1): 
       w1 = w0 * g1/g2 * Jm_R/Km_R  
    print(w1)

    fac = 100.0
    if (iunst == 0):
      ppert = fac*(g*rho/k)*w0*Jm * np.cos(k*Z + m*THETA - omega*t)
    elif (iunst == 1):
      ppert = fac*(g1*rho/k)*w0*Jm * np.cos(k*Z + m*THETA - omega*t)
  
    for i in np.arange(nmax):
      for j in np.arange(nmax):
        for kk in np.arange(nmax):

          if (RR[i,j, kk] > R):
            u[i,j,kk] = w1 *  Kpr[i,j, kk] *  np.sin(k*Z[i,j,kk] + m*THETA[i,j,kk] - omega*t)
            v[i,j,kk] = m*w1/(k*RR[i,j,kk]) * Km[i,j,kk] * np.cos(k*Z[i,j,kk] + m*THETA[i,j,kk] - omega*t)
            w[i,j,kk] = w1 * Km[i,j,kk] * np.cos(k*Z[i,j,kk] + m*THETA[i,j,kk] - omega*t)

            if (iunst == 0):  
              ppert[i,j,kk]  = fac*(g*rho/k)*w1*Km[i,j,kk] * np.cos(k*Z[i,j,kk] + m*THETA[i,j,kk] - omega*t)      
            elif (iunst == 1):
              ppert[i,j,kk]  = fac*(g2*rho/k)*w1*Km[i,j,kk] * np.cos(k*Z[i,j,kk] + m*THETA[i,j,kk] - omega*t)

            vv[i,j,kk]  = 0.8 * Omega* R**2/RR[i,j,kk] + v[i,j,kk]
            if (iretro == 1):
              vv [i,j,kk]  = 2.5 * Omega* R**2/RR[i,j,kk] + v[i,j,kk]

    # Outer  
  
    pbase = p0 - rho/2.0 * Omega**2 * R**4 / RR**2

    # Inner

    for i in np.arange(nmax):
      for j in np.arange(nmax):
        for kk in np.arange(nmax):
          if (RR[i,j, kk] <= R):
            pbase[i,j,kk] = p0 - rho/2.0 * Omega**2 * R**2 + rho/2.0 * Omega**2 * (RR[i,j,kk]**2 - R**2)

    p = pbase + ppert

    iso_val = np.min(np.min(p)) + 1.0
    verts, faces, dummy, dummy = measure.marching_cubes(p, iso_val, spacing=(0.1, 0.1, 0.1))


    xcoord = ((verts[:, 0] - 1.0) * 3.0/8.0 -xmax)
    ycoord = (verts[:,1] - 1.0) * 3.0/8.0 - ymax

    print('Plotting...')
    fig = plt.figure(figsize=(5,10), facecolor = 'w', edgecolor = 'k')
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=5, azim=45.0)
    ax.plot_trisurf(xcoord, ycoord, faces, verts[:, 2], antialiased=True) 
    ax.set_xticks([1.5, 0, -1.5])
    ax.set_yticks([1.5, 0, -1.5])
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_xlim(-1.5,1.5)
    ax.set_ylim(-1.5,1.5)
    ax.set_zlim(0,10)

    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel('z (m)', rotation=90)

    if (nt < 10):
      tstr = '00'+str(nt)
    elif (nt >= 10 and nt < 100):
      tstr = '0'+str(nt)
    elif (nt >= 100):
      tstr = str(nt)     


    if (iunst == 0):
      filename = "rankine_pure_m1_"+tstr+".png"
    elif (iunst == 1):
      filename = "rankine_axial_m1_"+tstr+".png"

    plt.savefig("./loops/"+filename) 
    print('Saved figure: ', filename)
    plt.show()
    plt.close(fig)

  quit()

#========================================================================
# SWIRLING RANKINE VORTEX (Loiseleux m = 1 wave)
#========================================================================

elif (icase == 3): # Rankine vortex

  nmax = 201
  R = 1.0
  xmax = R + 0.5
  ymax = R + 0.5
  zmin = 0.0
  zmax = 8.0
  ntmax = 100
  w0 = 0.22 # 1.0 # 0.22
  rho = 1.0
  p0 = 100000.0
  Omega = 1.0# 4.5 # 1.0
  W = 1.0

  m = 1
  k = 1.1
  omega = 1.8439990700321944 + 0.1449115402624638j

# -------------------
# Done with settings
# -------------------

  g1 = omega - m*Omega - W*k
  g2 = omega - m*Omega

  beta = k * np.sqrt(4.0*Omega**2/g1**2 - 1.0)

  for nt in np.arange(0, ntmax):

    print('Step: ', nt, ntmax)
    X,Y,Z = np.mgrid[(-nmax/2.0):(nmax/2.0):1, (-nmax/2.0):(nmax/2.0):1, 0:nmax:1]
    X = X * 2.0 * xmax/nmax
    Y = Y * 2.0 * ymax/nmax
    Z = Z * zmax/ nmax

    XX = X * np.cos(np.pi) + Y * np.sin(np.pi)
    YY = X * np.sin(np.pi) + Y * np.cos(np.pi)

    THETA =  np.pi + np.arctan2(YY,XX)

    RR = np.sqrt(X**2 + Y**2)

    Jm  = scipy.special.jv (m, beta* RR)
    Jpr = scipy.special.jvp(m, beta*RR)

    Km  = scipy.special.kv (m, k*RR)
    Kpr = scipy.special.kvp(m, k*RR)

    tmax = 2.0*pi/omega.real
    t   =  tmax/np.real(ntmax) * nt

# Normal mode

    argu = k*Z + m*THETA - omega*t
    wave = np.exp(argu*1j)

# Perturbation variables

    u = 1j * w0 * (k/beta * Jpr - k/(beta**2 * g1) * 2.0 * Omega * m / RR * Jm) * wave  
    v = w0 * k / (g1 * beta**2) * (2.0*Omega*beta*Jpr - (m*g1/RR)*Jm) * wave
    w = w0 * Jm * wave
    ppert = rho*g1/k * w

    pbase = p0 - 0.5*rho*Omega**2*R**2 + 0.5*rho*Omega**2*(RR**2 - R**2)
    vv = Omega * RR + v.real

    Jpr_R = scipy.special.jvp(m, beta*R)
    Jm_R  = scipy.special.jv (m, beta* R) 
    Kpr_R = scipy.special.kvp(m, k*R)
    Km_R  = scipy.special.kv (m, k*R)

  # Matching dispacement 

    w1 = - k * w0 *g2/(g1**2 * beta**2) * (g1*beta*Jpr_R - 2.0*m*Omega/R * Jm_R) / Kpr_R
    print(w1)

    w1 = w0 * g1/g2 * Jm_R/Km_R
    print(w1)

  # Outer region  

    for i in np.arange(nmax):
      for j in np.arange(nmax):
        for kk in np.arange(nmax):

          if (RR[i,j, kk] > R):
            u[i,j,kk] = -1j*w1 * Kpr[i,j,kk] * wave[i,j,kk]
            v[i,j,kk] = m*w1/(k*RR[i,j,kk]) * Km[i,j,kk] * wave[i,j,kk]
            w[i,j,kk] = w1 * Km[i,j,kk] * wave[i,j,kk]
            vv[i,j,kk]  = Omega* R**2/RR[i,j,kk] + v.real[i,j,kk]
            ppert[i,j,kk] = rho*g2/k * w[i,j,kk]
            pbase[i,j,kk] = p0 - 0.5*rho*Omega**2 * R**4/RR[i,j,kk]**2
  
    p = pbase + ppert.real
 
    print(np.min(ppert), np.min(p), np.max(ppert),np.max(p))

    iso_val = 99998.99563524182 # = initial np.min(p) + 0.015

    verts, faces, dummy, dummy = measure.marching_cubes(p, iso_val, spacing=(0.1, 0.1, 0.1))
    
    xcoord = ((verts[:, 0] - 8.0) * 3.0/4.0 -xmax)
    ycoord = (verts[:,1] - 8.0) * 3.0/4.0 - ymax
    zcoord = verts[:, 2]/20.0 * 8.0

    print('Plotting...')
    fig = plt.figure(figsize=(5,10), facecolor = 'w', edgecolor = 'k')
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=5, azim=45.0)
    ax.plot_trisurf(xcoord, ycoord, faces, zcoord, antialiased=True)
    ax.set_xticks([1.5, 0, -1.5])
    ax.set_yticks([1.5, 0, -1.5])
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_xlim(-1.5,1.5)
    ax.set_ylim(-1.5,1.5)
    ax.set_zlim(0,8)

    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel('z (m)', rotation=90)

    if (nt < 10):
      tstr = '00'+str(nt)
    elif (nt >= 10 and nt < 100):
      tstr = '0'+str(nt)
    elif (nt >= 100):
      tstr = str(nt)     

    filename = "rankine_axial_m1_"+tstr+".png"
    plt.savefig("./loops/"+filename) 
    print('Saved figure: ', filename)
   # plt.show()
    plt.close(fig)

  quit()

#=========================================================================
# HOLLOW VORTEX (Rotunno)
#=========================================================================

elif (icase == 4): # Rotunno's cylindrical vortex sheet

  m = 2
  nmax = 101
  R = 1.0
  W = 1.0
  xmax = R + 0.5
  ymax = R + 0.5
  zmin = 0.0
  zmax = 15.0
  p0 = 100000.0
  rho = 1.0
  ntmax = 100
  z = 0.0
  w0 = 1.0
  Omega = 1.0

  # From rotunno.py

  k = 1.0
  omega = 0.9086234646938594 + 1.6913247335984338j
  g1 = omega + W*k
  g2 = omega - m*Omega - W*k

  for nt in np.arange(0, ntmax):

    X,Y,Z = np.mgrid[(-nmax/2.0):(nmax/2.0):1, (-nmax/2.0):(nmax/2.0):1, 0:nmax:1]
    X = X * 2.0 * xmax/nmax
    Y = Y * 2.0 * ymax/nmax
    Z = Z * zmax/ nmax

    XX = X * np.cos(np.pi) + Y * np.sin(np.pi)
    YY = X * np.sin(np.pi) + Y * np.cos(np.pi)

    THETA =  np.pi + np.arctan2(YY,XX)

    RR = np.sqrt(X**2 + Y**2)

    Km  = scipy.special.kv (m, k*RR)
    Kpr = scipy.special.kvp(m, k*RR)

    Im  = scipy.special.iv (m, k*RR)
    Ipr = scipy.special.ivp(m, k*RR)

    tmax = 2.0*pi/omega.real
    t   =  tmax/np.real(ntmax) * nt

    argu = 1j*(k*Z + m*THETA - omega*t)
    wave = np.exp(argu)

    u_cp = -1j * w0 *  Ipr * wave
    v_cp = w0 * m/(RR * k) * Im * wave
    w_cp = w0 * Im *wave 

    u = u_cp.real
    v = v_cp.real
    w = w_cp.real

    ppert_cp = w0 * g1 * rho/k*Im * wave
    ppert = ppert_cp.real

    pbase_min = p0 - 0.5 * rho * Omega**2 * R**2 
    pbase = pbase_min + np.zeros(shape=(nmax,nmax,nmax))
    vv = np.zeros(shape=(nmax,nmax,nmax))

    Kpr_R = scipy.special.kvp(m, k*R)
    Ipr_R = scipy.special.ivp(m, k*R)

    w1 = w0 * Ipr_R/Kpr_R * g2/g1

    print(w1)

   # Overwrite values for the outer region

    for i in np.arange(nmax):
      for j in np.arange(nmax):
        for kk in np.arange(nmax):
          if (RR[i,j,kk] > R):
    
            u[i,j,kk] = -1j * w1 * Kpr[i,j,kk] * wave[i,j,kk]
            v[i,j,kk] = w1 * m/(k*RR[i,j,kk]) * Km[i,j,kk] * wave[i,j,kk]
            w[i,j,kk] = w1 * Km[i,j,kk] * wave[i,j,kk] 
            ppert[i,j,kk] = w1 * g2*rho/k * Km[i,j,kk]  * wave[i,j,kk]
            pbase[i,j,kk]  = p0 - rho/2.0 * Omega**2 * R**4/RR[i,j,kk]**2
            vv[i,j,kk]  = Omega* R**2/RR[i,j,kk]  # Just the base state 


    # Reducing the growth rate to 0.1 s-1 (about 6% of what it really is) or
    # the wave rapidly grows outside the domains boundary). 

    amplitude = np.exp(0.1*t) * np.exp((-omega.imag)*t)

    print('ampl:', amplitude, np.exp((omega.imag)*t), nt)

    p = pbase + amplitude * ppert.real

    print(amplitude*np.min(ppert), np.min(p), np.max(ppert),np.max(p))

    iso_val = 99999.3494215034 # np.min(np.min(p)) + 0.2

    verts, faces, dummy, dummy = measure.marching_cubes(p, iso_val, spacing=(0.1, 0.1, 0.1))

    xcoord = ((verts[:, 0] - 1.0) * 3.0/8.0 -xmax)
    ycoord = (verts[:,1] - 1.0) * 3.0/8.0 - ymax

    print('Plotting...')
    fig = plt.figure(figsize=(5,10), facecolor = 'w', edgecolor = 'k')
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=5, azim=45.0)
    ax.plot_trisurf(xcoord, ycoord, faces, verts[:, 2], antialiased=True)
    ax.set_xticks([1.5, 0, -1.5])
    ax.set_yticks([1.5, 0, -1.5])
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_xlim(-1.5,1.5)
    ax.set_ylim(-1.5,1.5)
    ax.set_zlim(0,10)

    ax.zaxis.set_rotate_label(False) 
    ax.set_zlabel('z (m)', rotation=90)

    if (nt < 10):
      tstr = '00'+str(nt)
    elif (nt >= 10 and nt < 100):
      tstr = '0'+str(nt)
    elif (nt >= 100):
      tstr = str(nt)

    filename = "rotunno_m2_"+tstr+".png"

    plt.savefig("./loops/"+filename)
    print('Saved figure: ', filename)
  #  plt.show()
    plt.close(fig)

  quit()

# All done.

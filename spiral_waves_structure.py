import numpy as np
import matplotlib.pyplot as plt 
import scipy
from scipy import optimize
from scipy import special

# Plot perturbation and total flows for spiral Kelvin modes
# in the horizontal plane.  Questions/comments: johannes.dahl@ttu.edu.

# Select case:

#icase = 1 # Bounded Couette flow
icase = 2 # Rankine vortex without axial flow
#icase = 3 #  Rankine vortex WITH axial core flow, stable modes
#icase = 4 # Rankine vortex WITH axial core flow, unstable modes
#icase = 5 # Cylindrical vortex sheet

#=========================================================================
# BOUNDED VORTEX (Couette)
#=========================================================================

if icase == 1: # Bounded
  nmax = 201
  R = 1.0
  xmax = R + 0.5 
  ymax = R + 0.5 
  k = 1.0 #2.0*np.pi# 1.0
  t = 0.0
  z = 0.0
  w0 = 1.0
  Omega = 1.0
  iretro = 0

# Bounded flow: Allowed eigenvalues

#  m = 0
#  alpha = 3.8317
# alpha = 7.0156

#  m = 1
#  alpha = 4.5590875
#  dr = 0.1

#  alpha = 7.7675125 
#  alpha = 10.935375 
#  dr = 0.2

# Retrograde modes

  iretro = 1
  m = 1
  k = 1.2
  alpha = 2.87207143
#  alpha = 6.09514286
#  alpha =9.26721429
  dr = 0.2

#  m = 2
#  alpha = 5.5794375 
#  alpha  = 8.86875
#   alpha = 12.074625
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

  x = np.linspace(-xmax,xmax,nmax) 
  y = np.linspace(-ymax,ymax,nmax)

  X,Y = np.meshgrid(x,y)

  XX = X * np.cos(np.pi) + Y * np.sin(np.pi)
  YY = X * np.sin(np.pi) + Y * np.cos(np.pi)

  THETA =  np.pi + np.arctan2(YY,XX)

  RR = np.sqrt(X**2 + Y**2)

# Calculate frequencies and velocities

  g = 2.0 * Omega / (np.sqrt(alpha**2 / (R*k)**2 + 1.0))

  if (iretro == 1):
    g = -2.0 * Omega / (np.sqrt(alpha**2 / (R*k)**2 + 1.0))
  
  omega = g + m*Omega
  print(g, omega)
 
  Jm  = scipy.special.jv (m, beta* RR)
  Jpr = scipy.special.jvp(m, beta*RR)

  u = -w0 * g * (g*beta*Jpr - 2.0*m*Omega/RR * Jm) / (k * (4.0 * Omega**2 - g**2)) * np.sin(k*z + m*THETA - omega*t)
  V = Omega * RR

  vv = Omega*RR + w0 * g * (2.0*Omega * beta * Jpr - m*g/RR * Jm) / (k * (4.0 * Omega**2 - g**2)) * np.cos(k*z + m*THETA - omega*t)
  v =  w0 * g * (2.0*Omega * beta * Jpr - m*g/RR * Jm) / (k * (4.0 * Omega**2 - g**2)) * np.cos(k*z + m*THETA - omega*t)

  for i in np.arange(nmax):
    for j in np.arange(nmax):
      if (RR[i,j] > R):
        u[i,j] = 0.0
        v[i,j] = 0.0
        vv[i,j] = 0.0

# Displacement of material loop

  Rm = R-dr

  Jpr_R = scipy.special.jvp(m, beta*Rm)
  Jm_R  = scipy.special.jv (m, beta* Rm)
  Kpr_R = scipy.special.kvp(m, k*Rm)
  Km_R  = scipy.special.kv (m, k*Rm)

  bd_scale = 0.5
  nn = 101
  theta = np.linspace(0,2.0*np.pi, nn)
  rpr = -bd_scale * w0/g * (k/beta * Jpr_R - 2.0 * Omega*m*k/(g*beta**2 * Rm) * Jm_R) * np.cos(k*z + m*theta - omega*t)

#  d = R + rpr 

  xpr = rpr*np.cos(theta)
  ypr = rpr*np.sin(theta)

  xr = Rm*np.cos(theta)
  yr = Rm*np.sin(theta)

  xd =  xr + xpr
  yd = yr + ypr

 # Cartesian velocity components (rotated oppositely to the normal sense because of
 # how Python defines the grid)

  uc = np.cos(THETA) * u - np.sin(THETA) * v
  vc = np.sin(THETA) * u + np.cos(THETA) * v

  uuc = np.cos(THETA) * u - np.sin(THETA) * vv
  vvc = np.sin(THETA) * u + np.cos(THETA) * vv

# ---------------------------------------------
# Plot stuff
# ---------------------------------------------

  plt.figure(figsize=(10,5), facecolor = 'w', edgecolor = 'k')

  xstart = 0
  xend = nmax-0
  ystart = 0
  yend = nmax-0
  istride = 4
  jstride = 4

  ax = plt.subplot(1,2,1)

  plt. contour(X,Y,RR, levels = [R], colors='k', linewidths=2)
  plt.plot(xd, yd, color='r', linewidth='2.0')

  arrowsize=0.6 
  plt.streamplot(X,Y,uc,vc, density=2.5,arrowsize=0.6, linewidth=0.3, color='k')
  plt.xlabel('x')
  plt.xlabel('y')

  # Right panel

  ax = plt.subplot(1,2,2)

  plt. contour(X,Y,RR, levels = [R], colors='k', linewidths=2)

  plt.quiver(X[xstart:xend:istride, ystart:yend:jstride],
  Y[xstart:xend:istride, ystart:yend:jstride],
  uuc[xstart:xend:istride, ystart:yend:jstride],
  vvc[xstart:xend:istride, ystart:yend:jstride], scale = 20.0)
  plt.xlabel('x')

  if (m == 1):
    plt.savefig("./vessel_horiz_structure_m1.eps", format = 'eps')
  elif (m == 4):
    plt.savefig("./vessel_horiz_structure_m4.eps", format = 'eps')

  plt.show()

#========================================================================
# PURE RANKINE VORTEX
#========================================================================

elif (icase == 2): # Rankine vortex

  nmax = 201
  R = 1.0
  xmax = R + 0.5
  ymax = R + 0.5
  t = 0.0
  z = 0.0
  w0 = 1.0  # 0.4 for retro if Omega = 1
  rho = 1.0
  Omega = 1.0

# Axisymmetric modes 
# -----------------
#
#  m = 0
#  z = np.pi # 0.0
#  k = 1.0
#  alpha = 2.650875 # Radial modes
#  alpha = 5.642
#  alpha = 8.73325

# ------------------------------
# Spiral modes (m = 1), cograde
# ------------------------------
#
#  m = 1
#  iretro = 0

#  k = 1.0
#  alpha = 3.957875 # fundamental
#  alpha = 7.09156 # second radial
#  alpha = 10.22756 # third radial

#  k = 2.0
#  alpha = 4.11936 # fundamental
#  alpha = 7.20132 # second radial
#  alpha = 10.3089 # third radial
#
# Double spiral modes, cograde
# ---------------------------

#  m = 2
#  k = 1.0
#  alpha = 5.2045
#  alpha = 8.4623
#  alpha =11.653425

# --------------------------------
# Spiral modes (m = 1), retrograde
# --------------------------------

# Countergrade (structureless), m = 1 
# -----------------------------------

#  m = 1
#  iretro = 1
#  k = 0.05
#  alpha = 0.086  # longwave limit displacement mode (not even shear vorticity)

#  k = 1.0
#  alpha = 1.13028 

#  k = 1.07
#  alpha= 1.17312 # structureless / displacement countergrade mode 

#  k = 1.1
#  alpha = 1.19064

# Retrograde modes, m = 1, k = 1.0
# --------------------------------
  m = 1
  iretro = 1
  
#  k = 1.0
#  alpha = 3.99984 # First retrograde mode (first sructured mode)
#  alpha = 7.10508 # Second radial mode
#  alpha = 10.2342 # Third radial mode

  k = 1.1
#  alpha = 4.0242 # First structured retrograde   
  alpha = 7.11852 # Second retrograde
#  alpha = 10.24332  # Third retrograde

# Retrograde double-helical modes (m = 2)
# -------------------------------------
#  iretro = 1
#  m = 2
#  k = 1.1
#  alpha = 1.58712  # structureless (countergrade)
#  alpha = 5.2488   # First strucgured retrograde
#  alpha = 8.48256   # Secodn structureless retrograde

# -------------------------------------------------------------------
# Done with settings
# -------------------------------------------------------------------

  beta = alpha/R

  x = np.linspace(-xmax,xmax,nmax)
  y = np.linspace(-ymax,ymax,nmax)

  X,Y = np.meshgrid(x,y)

  XX = X * np.cos(np.pi) + Y * np.sin(np.pi)
  YY = X * np.sin(np.pi) + Y * np.cos(np.pi)

  THETA =  np.pi + np.arctan2(YY,XX)

  RR = np.sqrt(X**2 + Y**2)

  # Intrinsic freq

  g = 2.0 * Omega / (np.sqrt(alpha**2 / (R*k)**2 + 1.0))

  # Modify for retrograde and unstable modes
 
  if (iretro == 1):
    g = -2.0 * Omega / (np.sqrt(alpha**2 / (R*k)**2 + 1.0))

  omega = g + m*Omega

  Jm  = scipy.special.jv (m, beta* RR)
  Jpr = scipy.special.jvp(m, beta*RR)

  Km  = scipy.special.kv (m, k*RR)
  Kpr = scipy.special.kvp(m, k*RR)

  u = w0 * g * (g*beta*Jpr - 2.0*m*Omega/RR * Jm)     / (k * (4.0 * Omega**2 - g**2)) * (-np.sin(k*z + m*THETA - omega*t))
  v = w0 * g * (2.0*Omega * beta * Jpr - m*g/RR * Jm) / (k * (4.0 * Omega**2 - g**2)) *   np.cos(k*z + m*THETA - omega*t)
  vv = 0.8*Omega * RR + v 
  ppert = w0 * g*rho/k * Jm * np.cos(k*z + m*THETA - omega*t) 

  if (iretro == 1):
    vv = 2.5*Omega * RR + v
  
  w = w0 * Jm * np.cos(k*z + m*THETA - omega*t)

  Jpr_R = scipy.special.jvp(m, beta*R)
  Jm_R  = scipy.special.jv (m, beta* R) 
  Kpr_R = scipy.special.kvp(m, k*R)
  Km_R  = scipy.special.kv (m, k*R)

  # Matching dispacement using kinematic BC

  w1 = - k * w0/(g * beta**2) * (g*beta*Jpr_R - 2.0*m*Omega/R * Jm_R) / Kpr_R
  print(w1)

  # calculate material boundary displacement
 
  nn = 101
  theta = np.linspace(0,2.0*np.pi, nn)

  rpr = w0/g * Kpr_R/Km_R * Jm_R * np.cos(k*z + m*theta - omega*t)

  xpr = rpr*np.cos(theta)
  ypr = rpr*np.sin(theta) 

  xr = R*np.cos(theta)
  yr = R*np.sin(theta)  

  dx =  xr + xpr
  dy = yr + ypr

  # Matching pressure (sanity check)

  w1 = w0 * Jm_R/Km_R
  
  print(w1)

  for i in np.arange(nmax):
    for j in np.arange(nmax):
      if (RR[i,j] > R):
        u[i,j] = w1 *  Kpr[i,j] *  np.sin(k*z + m*THETA[i,j] - omega*t)
        v[i,j] = m*w1/(k*RR[i,j]) * Km[i,j] * np.cos(k*z + m*THETA[i,j] - omega*t)
        w[i,j] = w1 * Km[i,j] * np.cos(k*z + m*THETA[i,j] - omega*t)       
        vv[i,j]  = 0.8 * Omega* R**2/RR[i,j] + v[i,j]
        ppert[i,j] = w1 * g * rho/k * Km[i,j] * np.cos(k*z + m*THETA[i,j] - omega*t)
        if (iretro == 1):
#          vv [i,j]  = Omega* R**2/RR[i,j] + v[i,j]     
          vv [i,j]  = 2.5 * Omega* R**2/RR[i,j] + v[i,j]

  uc = np.cos(THETA) * u - np.sin(THETA) * v
  vc = np.sin(THETA) * u + np.cos(THETA) * v

  uuc = np.cos(THETA) * u - np.sin(THETA) * vv
  vvc = np.sin(THETA) * u + np.cos(THETA) * vv

# -----------------
# Plot fields
# -----------------

  plt.figure(figsize=(10,5), facecolor = 'w', edgecolor = 'k')

  xstart = 0
  xend = nmax-0
  ystart = 0
  yend = nmax-0
  istride = 4
  jstride = 4

#----------------------------------------------
# Left panel
#----------------------------------------------

  ax = plt.subplot(1,2,1)
  
  plt.contour(X,Y,RR, levels = [R], colors='k', linestyles='dashed', linewidths=2)
  plt.streamplot(X,Y,uc,vc, density=3.0,arrowsize=0.6, linewidth=0.3, color='k')

  plt.xlim(-1.5, 1.5)
  plt.xlabel('x')
  plt.ylabel('y')

#----------------------------------------------
# Right panel
#----------------------------------------------

  ax = plt.subplot(1,2,2)

  plt.contour(X,Y,RR, levels = [R], colors='k', linestyles='dashed', linewidths=2)
  
  sc = 20.0
  if (iretro == 1):
    sc = 50.0
    istride = 6
    jstride = 6

  plt.quiver(X[xstart:xend:istride, ystart:yend:jstride],
  Y[xstart:xend:istride, ystart:yend:jstride],
  uuc[xstart:xend:istride, ystart:yend:jstride],
  vvc[xstart:xend:istride, ystart:yend:jstride], scale = sc)
  plt.xlabel('x')

  plt.savefig("./rankine_pure_horiz_structure_m1_third.eps", format = 'eps')

  plt.show()

#========================================================================
# LOISELEUX: Stable modes (Ranking vortex with swirling core)
#========================================================================

elif (icase == 3): # Rankine vortex

  nmax = 201
  rho = 1.0
  R = 1.0
  xmax = R + 0.5
  ymax = R + 0.5
  t = 0.0
  z = 0.0
  w0 = 1.0
  Omega = 1.0

# Stable retrograde spiral modes, m = 1

  m = 1
  W = 1.0
  
# k = 1.1
#--------

  k = 1.1 

# structureless:
#  gg = -1.8387840638293542  # These are all core intrinsic frequencies!
#  alpha = 0.4706117652941324

#  gg = -0.5215371946925622
#  alpha = 4.07235180879522 

  gg = -0.2646140266178536
  alpha = 8.240906022650567

# k = 2.2
#---------

#  k = 2.2

#  alpha = 0.19272981824545615  
#  gg = -1.9923693514076495

#  alpha = 4.8543213580339515 
#  gg = -0.8255808155862709
 
#  alpha = 8.58086452161304
#  gg = -0.49670373994524386

# -------------------
# Done with settings
# -------------------

  beta = alpha/R

  x = np.linspace(-xmax,xmax,nmax)
  y = np.linspace(-ymax,ymax,nmax)

  X,Y = np.meshgrid(x,y)

  XX = X * np.cos(np.pi) + Y * np.sin(np.pi)
  YY = X * np.sin(np.pi) + Y * np.cos(np.pi)

  THETA =  np.pi + np.arctan2(YY,XX)

  RR = np.sqrt(X**2 + Y**2)

  # Intrinsic frequency

  g1 = gg
  omega = g1 + m*Omega + W*k

  Jm  = scipy.special.jv (m, beta* RR)
  Jpr = scipy.special.jvp(m, beta*RR)

  Km  = scipy.special.kv (m, k*RR)
  Kpr = scipy.special.kvp(m, k*RR)

  u = w0 * g1 * (g1*beta*Jpr - 2.0*m*Omega/RR * Jm)     / (k * (4.0 * Omega**2 - g1**2)) * (-np.sin(k*z + m*THETA - omega*t))
  v = w0 * g1 * (2.0*Omega * beta * Jpr - m*g1/RR * Jm) / (k * (4.0 * Omega**2 - g1**2)) *   np.cos(k*z + m*THETA - omega*t)
  vv = Omega * RR + v 
  ppert = w0 * g1*rho/k * Jm * np.cos(k*z + m*THETA - omega*t) 

  vv = Omega * RR + v

  w = w0 * Jm * np.cos(k*z + m*THETA - omega*t)

  Jpr_R = scipy.special.jvp(m, beta*R)
  Jm_R  = scipy.special.jv (m, beta* R) 
  Kpr_R = scipy.special.kvp(m, k*R)
  Km_R  = scipy.special.kv (m, k*R)

  # From kinematic BC

  g2 = g1 + W*k  
  w1 = - k * w0 *g2/(g1**2 * beta**2) * (g1*beta*Jpr_R - 2.0*m*Omega/R * Jm_R) / Kpr_R
  print(w1)

  # From dynamic BC
  w1 = w0 * g1/g2 * Jm_R/Km_R
  print(w1)

  # calculate material boundary displacement
 
  nn = 101
  theta = np.linspace(0,2.0*np.pi, nn)

  rpr = w1/g2*Kpr_R * np.cos(k*z + m*theta - omega*t)

  xpr = rpr*np.cos(theta)
  ypr = rpr*np.sin(theta) 

  xr = R*np.cos(theta)
  yr = R*np.sin(theta)  

  dx =  xr + xpr
  dy = yr + ypr

 # Assign values vor vortex periphery 

  for i in np.arange(nmax):
    for j in np.arange(nmax):
      if (RR[i,j] > R):
        u[i,j] = w1 *  Kpr[i,j] *  np.sin(k*z + m*THETA[i,j] - omega*t)
        v[i,j] = m*w1/(k*RR[i,j]) * Km[i,j] * np.cos(k*z + m*THETA[i,j] - omega*t)
        w[i,j] = w1 * Km[i,j] * np.cos(k*z + m*THETA[i,j] - omega*t)       
        vv[i,j]  = Omega* R**2/RR[i,j] + v[i,j]
        ppert[i,j] = w1 * g2 * rho/k * Km[i,j] * np.cos(k*z + m*THETA[i,j] - omega*t)

  # Cartesian components

  uc = np.cos(THETA) * u - np.sin(THETA) * v
  vc = np.sin(THETA) * u + np.cos(THETA) * v

  uuc = np.cos(THETA) * u - np.sin(THETA) * vv
  vvc = np.sin(THETA) * u + np.cos(THETA) * vv

# -----------------
# Plot fields
# -----------------

  plt.figure(figsize=(10,5), facecolor = 'w', edgecolor = 'k')

  xstart = 0
  xend = nmax-0
  ystart = 0
  yend = nmax-0
  istride = 6
  jstride = 6

#------------------------------------
# Left panel
#------------------------------------

  ax = plt.subplot(1,2,1)
  
  plt.contour(X,Y,RR, levels = [R], colors='k', linestyles='dashed', linewidths=2)
  plt.streamplot(X,Y,uc,vc, density=3.0,arrowsize=0.6, linewidth=0.3, color='k')

  plt.xlim(-1.5, 1.5)
  plt.xlabel('x')
  plt.ylabel('y')

#------------------------------------
# Right panel
#------------------------------------

  ax = plt.subplot(1,2,2)

  plt.contour(X,Y,RR, levels = [R], colors='k', linestyles='dashed', linewidths=2)

  sc = 20
  plt.quiver(X[xstart:xend:istride, ystart:yend:jstride],
  Y[xstart:xend:istride, ystart:yend:jstride],
  uuc[xstart:xend:istride, ystart:yend:jstride],
  vvc[xstart:xend:istride, ystart:yend:jstride], scale = sc)
  plt.xlabel('x')

  if (m == 1):
    plt.savefig("./rankine_pure_horiz_structure_m1_fundamental.eps", format = 'eps')
  elif (m == 2):
    plt.savefig("./rankine_pure_horiz_structure_m2_fundamental.eps", format = 'eps')

  plt.show()

#=======================================================================================
# UNSTABLE MODES SWIRLING-CORE RANKINE VORTEX (LOISELEUX)
#=======================================================================================

elif (icase == 4): # Rankine vortex

  nmax = 201
  R = 1.0
  xmax = R + 0.5
  ymax = R + 0.5
  rho = 1.0
  p0 = 100000.0
  t = 0.0
  z = 0.0
  w0 = 0.22
  W = 1.0
  Omega = 1.0
  ampl = 1.0

  m = 1

# Unstable modes

  k = 1.1
  omega = 1.8439990700321944 + 0.1449115402624638j

# Unstable, k = 2.2

#  k = 2.16723 
#  omega = 2.389383328779293 + 0.5114843452090718j

# Unstable, k = 3 (no resonance here)
#  k = 3.0
#  omega = 2.691177684069123 + 0.9381911087689195j

# -------------------
# Done with settings
# -------------------

  x = np.linspace(-xmax,xmax,nmax)
  y = np.linspace(-ymax,ymax,nmax)

  X,Y = np.meshgrid(x,y)

  XX = X * np.cos(np.pi) + Y * np.sin(np.pi)
  YY = X * np.sin(np.pi) + Y * np.cos(np.pi)

  THETA =  np.pi + np.arctan2(YY,XX)

  RR = np.sqrt(X**2 + Y**2)

  # Intrinsic freq

  g1 = omega - m*Omega - W*k
  g2 = omega - m*Omega
  
  g1r = g1.real
  g1i = g1.imag

  g2r = g2.real
  g2i = g2.imag

  omegar = omega.real
  omegai = omega.imag

  beta = k * np.sqrt(4.0*Omega**2/g1**2 - 1.0)

  betar = beta.real
  betai = beta.imag

  Jm  = scipy.special.jv (m, beta* RR)
  Jpr = scipy.special.jvp(m, beta*RR)

  Km  = scipy.special.kv (m, k*RR)
  Kpr = scipy.special.kvp(m, k*RR)

# Use complex form and let Python do the work

  argu = k*z + m*THETA - omega*t
  wave = np.exp(argu*1j)

  u = 1j * w0 * (k/beta * Jpr - k/(beta**2 * g1) * 2.0 * Omega * m / RR * Jm) * wave  
  v = w0 * k / (g1 * beta**2) * (2.0*Omega*beta*Jpr - (m*g1/RR)*Jm) * wave
  vv = Omega * RR + v.real

  pbase = p0 - 0.5*rho*Omega**2*R**2 + 0.5*rho*Omega**2*(RR**2 - R**2)

  w = w0 * Jm * wave
  ppert = rho*g1/k * w

  Jpr_R = scipy.special.jvp(m, beta*R)
  Jm_R  = scipy.special.jv (m, beta* R) 
  Kpr_R = scipy.special.kvp(m, k*R)
  Km_R  = scipy.special.kv (m, k*R)

# Match solutions using kinematic BC

  w1 = - k * w0 *g2/(g1**2 * beta**2) * (g1*beta*Jpr_R - 2.0*m*Omega/R * Jm_R) / Kpr_R
  print(w1)

# Check if the result is the same using the dynamic BC

  w1 = w0 * g1/g2 * Jm_R/Km_R
  print(w1)

  # calculate material boundary displacement
  
  nn = 101
  theta = np.linspace(0,2.0*np.pi, nn)

  scalar_argu = 1j*(k*z + m*theta - omega*t)
  scalar_wave = np.exp(scalar_argu)
  rpr_cp = w1/g2*Kpr_R * scalar_wave 

  rpr = 0.2 *  rpr_cp.real

  xpr = rpr*np.cos(theta)
  ypr = rpr*np.sin(theta) 

  xr = R*np.cos(theta)
  yr = R*np.sin(theta)  

  dx =  xr + xpr
  dy = yr + ypr

# Outer region

  for i in np.arange(nmax):
    for j in np.arange(nmax):
      if (RR[i,j] > R):
        u[i,j] = -1j*w1 * Kpr[i,j] * wave[i,j]
        v[i,j] = m*w1/(k*RR[i,j]) * Km[i,j] * wave[i,j]
        w[i,j] = w1 * Km[i,j] * wave[i,j]
        vv[i,j]  = Omega* R**2/RR[i,j] + v.real[i,j]
        ppert[i,j] = rho*g2/k * w[i,j]
        pbase[i,j] = p0 - 0.5*rho*Omega**2 * R**4/RR[i,j]**2

  uc = np.cos(THETA) * u.real - np.sin(THETA) * v.real
  vc = np.sin(THETA) * u.real + np.cos(THETA) * v.real

  p = pbase + ppert.real

  uuc = np.cos(THETA) * u.real - np.sin(THETA) * vv
  vvc = np.sin(THETA) * u.real + np.cos(THETA) * vv

# -----------------
# Plot fields
# -----------------

  plt.figure(figsize=(13,5), facecolor = 'w', edgecolor = 'k')

  xstart = 0
  xend = nmax-0
  ystart = 0
  yend = nmax-0
  istride = 6
  jstride = 6

  # -------------------------------
  # Left panel
  # -------------------------------

  ax = plt.subplot(1,2,1)
  
  plt.contour(X,Y,RR, levels = [R], colors='k', linestyles='dashed', linewidths=2)

  plt.streamplot(X,Y,uc,vc, density=3.0,arrowsize=0.6, linewidth=0.3, color='k')
  levs = np.linspace(-0.4, 0.4, 17)
  plt.contourf(X,Y,ppert, levels=levs, cmap='BrBG_r')
  plt.colorbar()

  plt.xlim(-1.5, 1.5)
  plt.xlabel('x')
  plt.ylabel('y')

  # -------------------------------
  # Right panel
  # -------------------------------

  ax = plt.subplot(1,2,2)

  plt.contour(X,Y,RR, levels = [R], colors='k', linestyles='dashed', linewidths=2)
  
  levs = np.linspace(-1.1,0, 23)
  plt.contourf(X,Y, p-p0, levels=levs, cmap='Blues_r')
  plt.colorbar() 

  sc = 15.0 

  plt.quiver(X[xstart:xend:istride, ystart:yend:jstride],
  Y[xstart:xend:istride, ystart:yend:jstride],
  uuc[xstart:xend:istride, ystart:yend:jstride],
  vvc[xstart:xend:istride, ystart:yend:jstride], scale = sc)
  plt.xlabel('x')

  plt.savefig("./rankine_axial_unstable_structure_m1.eps", format = 'eps')

  plt.show()
  quit()

#=========================================================================
# HOLLOW VORTEX (Rotunno)
#=========================================================================

elif (icase == 5): # Rotunno's cylindrical vortex sheet
 
  W = 1.0
  m = 2
  k = 1.0
  nmax = 201
  R = 1.0
  xmax = R + 0.5
  ymax = R + 0.5
  t = 0.0
  rho = 1.0
  z = 0.0
  w0 = 1.0
  Omega = 1.0
  dr = 0.0
  p0 = 100000.0
  # From rotunno.py

  omega = 0.9086234646938594 + 1.6913247335984338j

  omegar = omega.real
  omegai = omega.imag
  g1 = omega + W*k  # inner (V = Omega = 0; W < 0)
  g2 = omega - m*Omega - W*k

  g1r = g1.real
  g2r =  g2.real

  g1i = omegai
  g2i = omegai

  x = np.linspace(-xmax,xmax,nmax)
  y = np.linspace(-ymax,ymax,nmax)

  X,Y = np.meshgrid(x,y)

  XX = X * np.cos(np.pi) + Y * np.sin(np.pi)
  YY = X * np.sin(np.pi) + Y * np.cos(np.pi)

  THETA =  np.pi + np.arctan2(YY,XX)

  RR = np.sqrt(X**2 + Y**2)

  Km  = scipy.special.kv (m, k*RR)
  Kpr = scipy.special.kvp(m, k*RR)

  Im  = scipy.special.iv (m, k*RR)
  Ipr = scipy.special.ivp(m, k*RR)

  argu = 1j*(k*z + m*THETA - omega*t)
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
  pbase = pbase_min + np.zeros(shape=(nmax,nmax))

  vv = np.zeros(shape=(nmax,nmax))
  uu = np.zeros(shape=(nmax,nmax))

  # Boundary displacement

  Kpr_R = scipy.special.kvp(m, k*R)
  Ipr_R = scipy.special.ivp(m, k*R)

  Im_R = scipy.special.iv(m, k*R)
  Km_R = scipy.special.kv(m, k*R)

  # Should be identical for given k, omega
  
  w1 = w0 * Ipr_R/Kpr_R * g2/g1  # Kinematic
  print(w1, k)
  w1 = w0 * (g1 * Im_R - k*Omega**2*R/g1 * Ipr_R)/(g2 * Km_R)  # Dynamic
  print(w1, k)

  w1r = w1.real
  w1i = w1.imag

  C = w1r*g2r + w1i*g2i
  D = w1i*g2r - w1r*g2i

  nn = 101
  theta = np.linspace(0,2.0*np.pi, nn)
  rpr =  w0 /(g1r**2 + g1i**2) * Ipr_R * (g1r * np.cos(k*z + m*theta - omegar*t) + g1i*np.sin(k*z + m*theta - omegar*t))

#  d = R + rpr

  xpr = rpr*np.cos(theta)
  ypr = rpr*np.sin(theta)

  xr = R*np.cos(theta)
  yr = R*np.sin(theta)

  dx =  xr + xpr
  dy = yr + ypr

  A = w1r*g2r - w1i*g2i
  B = w1r*g2i + w1i*g2r

  for i in np.arange(nmax):
    for j in np.arange(nmax):
      if (RR[i,j] > R):

        u[i,j] = -1j * w1 * Kpr[i,j] * wave[i,j]
        v[i,j] = w1 * m/(k*RR[i,j]) * Km[i,j] * wave[i,j]
        w[i,j] = w1 * Km[i,j] * wave[i,j] 
        ppert[i,j] = w1 * g2*rho/k * Km[i,j]  * wave[i,j]
        pbase[i,j]  = p0 - rho/2.0 * Omega**2 * R**4/RR[i,j]**2
        vv[i,j]  = Omega* R**2/RR[i,j]  # Just the base state 
 
  p = pbase + ppert.real
 
  plt.plot(YY[:,101], ppert[:,101])
  plt.show()
  quit()

 # Cartesian components

  uc = np.cos(THETA) * u.real - np.sin(THETA) * v.real
  vc = np.sin(THETA) * u.real + np.cos(THETA) * v.real

# Base state

  uuc = np.cos(THETA) * uu - np.sin(THETA) * vv
  vvc = np.sin(THETA) * uu + np.cos(THETA) * vv

# Total

  utc = uuc + uc
  vtc = vvc + vc

# -----------------
# Plot fields
# -----------------

  plt.figure(figsize=(13,5), facecolor = 'w', edgecolor = 'k')

  xstart = 0
  xend = nmax-0
  ystart = 0
  yend = nmax-0
  istride = 6
  jstride = 6

  # --------------------------------------------------  
  # Left panel
  # --------------------------------------------------

  ax = plt.subplot(1,2,1)

  plt. contour(X,Y,RR, levels = [R], linestyles='dashed', colors='k', linewidths=2)
  
  levs = np.linspace(-0.4, 0.4, 17)
  cs=plt.contourf(X,Y,ppert.real, levels=levs, cmap='BrBG_r')
  plt.colorbar(cs) #, orientation='horizontal')

  arrowsize=0.6 
  plt.streamplot(X,Y,uc,vc, density=2.0,arrowsize=0.6, linewidth=0.3, color='k')
  plt.xlabel('x (m)')
  plt.ylabel('y (m)')

  # --------------------------------------------------  
  # Right panel
  # --------------------------------------------------

  ax = plt.subplot(1,2,2)

  plt.contour(X,Y,RR, levels = [R], colors='k', linestyles='dashed', linewidths=2)
  plt.plot(dx,dy,linewidth='2', color='r')

  levs = np.linspace(-0.9, 0, 19)
  cs=plt.contourf(X,Y,p-p0, cmap='Blues_r', levels = levs)
  plt.colorbar(cs)

  plt.quiver(X[xstart:xend:istride, ystart:yend:jstride],
  Y[xstart:xend:istride, ystart:yend:jstride],
  utc[xstart:xend:istride, ystart:yend:jstride],
  vtc[xstart:xend:istride, ystart:yend:jstride], scale = 10.0)

  plt.xlabel('x (m)')
  
  plt.savefig("./rotunno_m2.eps", format = 'eps')
  plt.show()

# -------------------------------------------------- 
# The End.
# --------------------------------------------------  




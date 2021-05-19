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

# Run in background: nohup python -u kelvin_dispersion.py &

#------------------------------------------------
# Set parameters
#------------------------------------------------

itest = 0

m = 2
Omega = 1.0
W =  1.0
nmax = 301
nnk = 31
kmax = 5.0
omega_max = 10.0

# Only for write-out during program execution
eps = 0.001 #  0.005 # 1.E-3 This value will depend on the resolution (nmax)

nnk = 6

#-------------------------------------------------
# Start
#-------------------------------------------------

itest = 0

if (itest == 1):
 
  m = 1
  W = 1.0
  dWdr =  1.E12 
  Omega = 1.0
  omega = 1.0

  k = 1.0

  nmax = 10
  omega_max = 5.0

  for mm in np.arange(nmax):  # incrementing real argument [0,2]
    for nn in np.arange(nmax):  # incrementing imaginary argument [0,2]
      omega_real = omega_max/(np.real(nmax-1)) * np.real(mm)
      omega_imag = omega_max/(np.real(nmax-1)) * np.real(nn) * 1j
      omega = omega_real + omega_imag

      g1 = omega - k*W
      g2 = omega - m * Omega + k*W 
  
      I_m = scipy.special.iv(m,k)
      K_m = scipy.special.kv(m,k)

      I_mp = scipy.special.ivp(m, k)
      K_mp = scipy.special.kvp(m, k)

      lhs = (g2/g1)**2 * I_mp/I_m - K_mp/K_m
      rhs = ((g1-g2)/g1**2) * I_mp/I_m * K_mp/K_m * dWdr

      print()
      print('STEP: ', mm, nn, omega, k)
      print('g', g1, g2, (g1-g2)/g1**2, g2**2/g1**2)
      print('Bessel', I_mp/I_mp, K_mp/K_m, I_mp/I_mp * K_mp/K_m)
      print('LHS RHS', lhs, rhs) 

  lhs = (g2/g1)**2 * I_mp/I_m - K_mp/K_m
  rhs = ((g1-g2)/g1**2) * I_mp/I_m * K_mp/K_m * dWdr


  plt.plot(k,lhs)
  plt.plot(k,rhs)

  plt.show()
  quit()

elif (itest == 2):

  m = 3
  S = 1.0 #  1.0
  R = 1.0
  k = np.linspace(0,10,101)
  eps = 1.0E-12
  alpha = scipy.special.iv(m,k*R) / (1.E-12 + k*R*scipy.special.iv(m-1,k*R) - m*scipy.special.iv(m,k*R)) 
  beta  = scipy.special.kv(m,k*R) / (1.E-12 + k*R*scipy.special.kv(m-1,k*R) + m*scipy.special.kv(m,k*R))
  t1 = -1.j/ (alpha+beta+eps) * (k * (beta-alpha) * S + beta*m)
  t2 = (k**2 * S**2 + 2.0*beta/(alpha + beta+eps) * m * k * S + beta/(alpha+beta+eps) * m**2)
  t3 = -1.0/(alpha+beta+eps) * 0.0 - 1.0/(alpha+beta+eps)**2 * (k*(beta-alpha)*S + beta*m)**2
  sarr  = np.sqrt(t2+t3)
  s  = np.zeros(101)

  print(t2)  

  for kk in np.arange(101):
    if ((t2[kk]+t3[kk]) >=0):

      alpha = scipy.special.iv(m,kk*R) / (kk*R*scipy.special.iv(m-1,kk*R) - m*scipy.special.iv(m,kk*R)+eps)
      beta  = scipy.special.kv(m,kk*R) / (kk*R*scipy.special.kv(m-1,kk*R) + m*scipy.special.kv(m,kk*R)+eps)
      t2[kk] = kk**2 * S**2 + 2.0*beta/(alpha + beta+eps) * m * kk * S + beta/(alpha+beta+eps) * m**2
 #   t3[kk] = -1.0/(alpha+beta) - 1.0/(alpha+beta)**2 * (kk*(beta-alpha)*S + beta*m)**2
      t3[kk] = -1.0/(alpha+beta+eps) * 0 - 1.0/(alpha+beta+eps)**2 * (kk*(beta-alpha)*S + beta*m)**2

      s[kk] = np.sqrt(t2[kk] + t3[kk])

  plt.plot(k,sarr)

  plt.xlim(0,10)
  plt.ylim(0,10)
  plt.show()
  quit()


  x = np.linspace(0,20,101)

  X,Y = np.meshgrid(x, x)

  Z  = X**2 + Y**2
  W = X**2 - Y**2
 
  plt.figure()

  C1 = plt.contour(X, Y, Z, [100], colors='r')
  C2 = plt.contour(X, Y, W, [0], colors='k')

  plt.figure()

#  A = [[1, 2], [4, 5]]
#  i = 1
#  j = 1
#  print(A, i,j, A[i][j]) # First index: which block; 2nd index: Element in block;
                         # Start counting at zero!

  # Get arrays containing the coordinates of contour in first plot
 
  aa = C1.allsegs[0]
  n1 = len(aa[0])
  arr1 = np.zeros(shape=(n1,2))
  
  for i in np.arange(n1):
    ar = aa[0][i]
    arr1[i,0] = ar[0]
    arr1[i,1] = ar[1]
  
  # Get arrays containing the coordinates of contour in second plot

  bb = C2.allsegs[0]
  n2 = len(bb[0])
  arr2 = np.zeros(shape=(n2,2))

  for j in np.arange(n2):
    ar = bb[0][j]
    arr2[j,0] = ar[0]
    arr2[j,1] = ar[1]

#  arr1[n,0], arr1[n,1] contains line segment starting and end points (of contour 1)

  for i in np.arange(n1-1): # Go through all segments of first line
    A = Point(arr1[i,0], arr1[i,1])
    B = Point(arr1[i+1,0], arr1[i+1,1])
    segm1 = LineString([A,B])
  
    for j in np.arange(n2-1):  # Go through all segments of second line
      C = Point(arr2[j,0], arr2  [j,  1])
      D = Point(arr2[j+1,0], arr2[j+1,1])
      segm2 = LineString([C,D])

      # Find intersecting coordinates (all except one being empty)

      int_pt = segm1.intersection(segm2)

      # Find the nonzero intersection points

      if not int_pt.is_empty:
        px = int_pt.x
        py = int_pt.y
        print(px,py)

  for ii, seg in enumerate(C2.allsegs[0]):
      plt.plot(seg[:,0], seg[:,1], '.-', label=ii)
  plt.legend(fontsize=9, loc='best')

  for ii, seg in enumerate(C1.allsegs[0]):
      plt.plot(seg[:,0], seg[:,1], '.-', label=ii)
  plt.legend(fontsize=9, loc='best')

  plt.show()
  quit()

#----------------------------------------------------------
# End test
#----------------------------------------------------------

nk = 0

or_arr = np.zeros(nnk)
oi_arr = np.zeros(nnk)
k_arr  = np.zeros(nnk)

for k in np.linspace(0,kmax,nnk):
 
  if (k == 0):
    k = 1.0E-12 # avoid k = 0

  # Initialize values of interest for complex omega

  fr2d = np.zeros(shape=(nmax,nmax))
  fi2d = np.zeros(shape=(nmax,nmax))

  counter = 0

  for mm in np.arange(nmax):  # incrementing real argument [0,2]
    for nn in np.arange(nmax):  # incrementing imaginary argument [0,2]
      omega_real = omega_max/(np.real(nmax-1)) * np.real(mm) 
      omega_imag = omega_max/(np.real(nmax-1)) * np.real(nn) * 1j
      omega = omega_real + omega_imag
   
# Tell us how far along we are

      if (nn == 50 and (mm % 50 == 0)):
        print(k, nn, mm, omega)

# Bessel function (1st kind) and its derivative

      I_m     = scipy.special.iv (m, k)
      I_mp    = scipy.special.ivp(m, k)

# Modified Bessel function (2nd kind) and its derivative

      K_m     = scipy.special.kv (m, k)
      K_mp    = scipy.special.kvp(m, k)
  
# Dispersion relation

      g1 = omega + W*k + 1.E-12   # inner (V = Omega = 0; W < 0)
      g2 = omega - m*Omega - W*k  # outer (W > 0)

      R = 1.0

      lhs = Omega**2 * R/(g1*g2) * I_mp * K_mp 
      rhs =-g2/(g1*k) * I_mp * K_m + g1/(g2*k) * K_mp * I_m

      D = lhs - rhs

      fr2d[mm,nn] = D.real
      fi2d[mm,nn] = D.imag

#  fig, ax = plt.subplots()

  x = np.linspace(0,omega_max,nmax)
  X,Y = np.meshgrid(x,x)

  fr2d = fr2d.T
  fi2d = fi2d.T

  #clevels= np.linspace(-2.0, 2.0, 20)
  #CS1 = plt.contour(X, Y, fr2d, levels = clevels) #  colors = 'k')
  #ax.clabel(CS1, inline=1, fontsize=10, fmt = '%1.1f')
  C1 = plt.contour(X, Y, fr2d, levels = [0.0], linewidths = 2.0, colors = 'k')

  #clevels= np.linspace(-2.0, 2.0, 20)
  #CS1 = plt.contour(X, Y, fi2d)#, levels = clevels) #  colors = 'k')
  #ax.clabel(CS1, inline=1, fontsize=10, fmt = '%1.1f')
  C2 = plt.contour(X, Y, fi2d, levels = [0.0], linewidths = 2.0, colors = 'r')

  plt.savefig("./dispersion_rotunno_2d.png", dpi=120)

#---------------------------------------
# Find intersection of zero-contours
#---------------------------------------

# Nested loops: Go through lines, then line segments of each plot

  aa = C1.allsegs[0]
  bb = C2.allsegs[0]

# How many contours there are (some insersecting with the axes)

  ncts1 = len(aa) # number of "blocks" or pairs defining a contour in first plot
  ncts2 = len(bb)

# See the structure
#print(aa)
#print(aa[0][0], aa[0][1], aa[0][2]) # First block of pairs
#print(aa[1][0], aa[1][1], aa[1][2]) # second bloc of pairs

  cc = 0
  px = np.zeros(100)
  py = np.zeros(100)

  for nc1 in np.arange(ncts1):

  # Number of data pairs of each line

    n1 = len(aa[nc1][:])

    arr1 = np.zeros(shape=(n1,2))

    for i in np.arange(n1):
      ar = aa[nc1][i]
      arr1[i,0] = ar[0]
      arr1[i,1] = ar[1]
  
  # Get arrays containing the coordinates of contour in second plot

    for nc2 in np.arange(ncts2):

  # Number of data pairs of each line

      n2 = len(bb[nc2][:])

      arr2 = np.zeros(shape=(n2,2))

      for j in np.arange(n2):
        ar = bb[nc2][j]
        arr2[j,0] = ar[0]
        arr2[j,1] = ar[1]

#  arr1[n,0], arr1[n,1] contains line segment starting and end points (of contour 1)

      for i in np.arange(n1-1): # Go through all segments of first line
        A = Point(arr1[i,0], arr1[i,1])
        B = Point(arr1[i+1,0], arr1[i+1,1])
        segm1 = LineString([A,B])
 
        for j in np.arange(n2-1):  # Go through all segments of second line
          C = Point(arr2[j,0], arr2  [j,  1])
          D = Point(arr2[j+1,0], arr2[j+1,1])
          segm2 = LineString([C,D])

        # Find intersecting coordinates (all except one being empty)

          int_pt = segm1.intersection(segm2)

    # Find the nonzero intersection points

          if (not int_pt.is_empty):
            px[cc] = int_pt.x
            py[cc] = int_pt.y + 1.0E-12 # So that the max function works
            ppx = int_pt.x
            ppy = int_pt.y

            cc = cc + 1
             
            # Update the arrays, always using the max value
             
            maxindex = np.argmax(py)
 
            or_arr[nk] = px[maxindex]
            oi_arr[nk] = py[maxindex]

            if (ppy > eps): # Only pick the unstable mode (omega_i > 0) 
              or_arr[nk] = ppx
              oi_arr[nk] = ppy
              print('')
              print('W ', 'step ', 'k ', 'omega_r ', 'omega_i ', W, nk, k, ppx, ppy)
              print(px[maxindex],py[maxindex])

              print(ppx + W*k + 1.E-12)   # inner (V = Omega = 0; W < 0)
              print(ppx - m*Omega - W*k) 

##^   
  nk = nk + 1

#------------------------------------------------------------
# Plot the dispersion relationship
#------------------------------------------------------------

# Remove Doppler effect

k_arr = np.linspace(0,kmax,nnk)

#or_arr = or_arr - m*S - (1.0 + a) * k_arr
#or_arr = or_arr - (m*S + (1.0 + a)*k_arr)

for i in np.arange(nnk):
  print(k_arr[i], or_arr[i], oi_arr[i])

if (1==0):
  R = 1.0
  S = 1
  k = np.linspace(0,kmax,nnk)
  eps = 1.0E-12
  alpha = scipy.special.iv(m,k*R) / (k*R*scipy.special.iv(m-1,k*R) - m*scipy.special.iv(m,k*R))
  beta  = scipy.special.kv(m,k*R) / (k*R*scipy.special.kv(m-1,k*R) + m*scipy.special.kv(m,k*R))
  t1 = -1.j/ (alpha+beta+eps) * (k * (beta-alpha) * S + beta*m)
  t2 = (k**2 * S**2 + 2.0*beta/(alpha + beta+eps) * m * k * S + beta/(alpha+beta+eps) * m**2)
  t3m = - 1.0/(alpha+beta+eps)**2 * (k*(beta-alpha)*S + beta*m)**2
  t3 =   -1.0/(alpha+beta+eps) - 1.0/(alpha+beta+eps)**2 * (k*(beta-alpha)*S + beta*m)**2
  

  sarr1 = np.sqrt(t2+t3m)
  sarr = np.sqrt(t2+t3)



plt.clf() # Somehow Python wants to plot the integrated previous figure

#plt.plot(k_arr, or_arr, 'k', linewidth = 2.0)
plt.plot(k_arr, oi_arr, 'r', linewidth = 2.0)
#plt.plot(k_arr, sarr1, 'b', linewidth = 2.0)
#plt.plot(k_arr, sarr, 'g', linewidth = 2.0)

plt.xlim(0, kmax)
plt.ylim(0,5.0)
plt.savefig("./dispersion.png", dpi=120)
plt.show()


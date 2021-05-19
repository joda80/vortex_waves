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

# Solving Loiseloix nondimensional dispersion equation
# Here we prescribe complex omega values and determine
# the complex beta; this is inserted into the dispersion
# relation and its roots are found in the complex plane.
# This will run for several hours.
# Run in background: nohup python -u kelvin_dispersion.py &

# Caveat: For large k, omega it doesn't work because
# the omega increases (Doppler effect due to nonzero W)
# but when Doppler-shifting the frequency they make no
# sense either.  Thought: Shifting it wong. Not noticing
# this effect for small'ish k, because the effect of
# the Doppler shift is small.

#------------------------------------------------
# Set parameters
#------------------------------------------------

S = 1.0
a = 0.0

# Production plot setting m = 0
#m = 0
#kmax = 2.0
#omega_max = 2.0
#nmax = 1501
#itest = 0
#nnk = 21 # For how many k

#m = +1
#------

S = 1.0
a = 0.0
m = 1
nmax = 601
nnk = 31
kmin = 0.001
kmax = 5.0

# Playing around

nnk = 5
kmin = 0.001
kmax = 0.25


# To find beta
nmax = 601
nnk = 2
kmax = 1.1 #1.0

omega_min_r = 0.0
omega_max_r = 4.0 # Use half of the value desired
omega_min_i = 0.0
omega_max_i = 8.0

# m = -1 
#-------

#S = 1.0
#a = 0.0
#m = -1
#nmax = 601
#nnk = 31
#kmin = 0.001
#kmax = 5.0

#For 2D plot

#nnk = 6
#nmax = 601

#omega_min_r = -5.0
#omega_max_r = 5.0
#omega_min_i = 0.0
#omega_max_i = 5.0

# Only for write-out during program execution
eps = 0.001 # 0.001# 0.001 #  0.005 # 1.E-3 This value will depend on the resolution (nmax)

#-------------------------------------------------
# Start
#-------------------------------------------------

nk = 0

or_arr = np.zeros(nnk)
oi_arr = np.zeros(nnk)
k_arr  = np.zeros(nnk)


for k in np.linspace(kmin,kmax,nnk):


  print('STEP', nk, k) 
  if (k == 0):
    k = 1.0E-12 # avoid k = 0

  # Initialize values of interest for complex omega

  fr2d = np.zeros(shape=(nmax,nmax))
  fi2d = np.zeros(shape=(nmax,nmax))
  X = np.zeros(shape=(nmax,nmax))
  Y = np.zeros(shape=(nmax,nmax))

  counter = 0

  for mm in np.arange(nmax):  # incrementing real argument [0,2]
    for nn in np.arange(nmax):  # incrementing imaginary argument [0,2]
      omega_real =  omega_min_r + 2.0*omega_max_r/(np.real(nmax-1)) * np.real(mm) 
      omega_imag = (omega_min_i + omega_max_i/(np.real(nmax-1)) * np.real(nn) ) * 1j
      omega = omega_real + omega_imag
  
      X[nn,mm] = omega_min_r + 2.0*omega_max_r/(np.real(nmax-1)) * np.real(mm)
      Y[nn,mm] = omega_min_i + omega_max_i/(np.real(nmax-1)) * np.real(nn)
 
# Tell us how far along we are

#      if (nn == 50 and (mm % 50 == 0)):
#        print(k, nn, mm, omega)

# Dimensionless coefficients:

      g    = omega - m*S - k
      beta = k * np.sqrt( (4.0*(S**2/g**2) - 1.0 ) )

# Bessel function (1st kind) and its derivative

      bess     = scipy.special.jv (m, beta)
      bess_der = scipy.special.jvp(m, beta)

# Modified Bessel function (2nd kind) and its derivative

      mbess     = scipy.special.kv (m, k)
      mbess_der = scipy.special.kvp(m, k)

# Dispersion relation

      lhs = (g+k)**2 * ( beta * bess_der/(bess+1.0E-12) - 2.0*S*m/g )
      rhs = -g**2 * beta**2 / k * mbess_der/(mbess+1.0E-12) 

      D = lhs - rhs

      fr2d[mm,nn] = D.real
      fi2d[mm,nn] = D.imag

#  fig, ax = plt.subplots()

  print('Done with omega-loop')    

#  x = np.linspace(omega_min_r, omega_max,nmax)
#  y = np.linspace(omega_min_i, omega_max,nmax)
#
#  X,Y = np.meshgrid(x,y)

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

  plt.savefig("./dispersion_2d.png", dpi=120)

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
#              or_arr[nk] = ppx
#              oi_arr[nk] = ppy
              print('')
              print('S ', 'step ', 'm', 'k ', 'omega_r ', 'omega_i ', S, nk, m, k, ppx, ppy)
              print(px[maxindex],py[maxindex])

##^   
  nk = nk + 1

#------------------------------------------------------------
# Plot the dispersion relationship
#------------------------------------------------------------

k_arr = np.linspace(0.001,kmax,nnk)

#or_arr = or_arr - m*S - (1.0 + a) * k_arr
#or_arr = or_arr - (m*S + (1.0 + a)*k_arr)

beta_arr = np.zeros(nnk)

for i in np.arange(nnk):

  g    = or_arr[i] - m*S - k_arr[i]
  beta_arr[i] = k_arr[i] * np.sqrt( (4.0*(S**2/g**2) - 1.0 ) )  

  print(g) 
  print('Output') 
  print(k_arr[i], or_arr[i], oi_arr[i], beta_arr[i])

plt.clf() # Somehow Python wants to plot the integrated previous figure

plt.plot(k_arr, or_arr, 'k', linewidth = 2.0)
plt.plot(k_arr, oi_arr, 'r', linewidth = 2.0)
plt.savefig("./dispersion.png", dpi=120)
plt.show()


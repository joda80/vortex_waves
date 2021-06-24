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

# Solving the nondimensional dispersion relation by 
# Loiseleux et al. (1998, PoF).
# The complex omega values are prescribed and determine
# the complex beta; this is inserted into the dispersion
# relation and its roots are found graphically in the complex plane.
#
# When no unstable growth rate is found, this is often tied
# to the resolution of the D field, so increasing nmax usually
# helps; also, omega_min and omega_max must cover the relevant domain
# (omega_min may be negative in some cases).  The present settings
# work, but nmax, omega_min, and omega_max may need to be adjusted
# if different azimuthal wavenumbers or axial wavenumbers are considered.
#
# Contact johannes.dahl@ttu.edu for questions.

#------------------------------------------------
# Set parameters
#------------------------------------------------

S = 1.0  # Swirl parameter
a = 0.0  # Co-flow parameter 

# m = 0
# -----

#m = 0 # Reproducing Loiseleux' plot (not used in teh MWR review paper)
#kmax = 2.0
#omega_max = 2.0
#nmax = 1501
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

omega_min_r = 0.0
omega_max_r = 4.0
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

#omega_min_r = -5.0
#omega_max_r = 5.0
#omega_min_i = 0.0
#omega_max_i = 5.0

# To make sure the solver only picks unstable modes, omega_i must be larger than eps

eps = 0.001 # This value needs to be larger when nmax is reduced

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

  counter = 0

  omega_r = np.linspace(omega_min_r, omega_max_r, nmax)
  omega_i = np.linspace(omega_min_i, omega_max_i, nmax) * 1j

  X, Y = np.meshgrid(omega_r, omega_i)

  Omega = X + Y  # Not to be confused with the base state angular velocity, which does
                 # not appear in the dimensionless formulation

  g    = Omega - m*S - k
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

  fr2d = D.real
  fi2d = D.imag

# Graphical solution (do not delete this part; C1 and C2 are used below)

  ox = np.linspace(omega_min_r, omega_max_r, nmax)
  oy = np.linspace(omega_min_i, omega_max_i, nmax)

  XX, YY = np.meshgrid(ox,oy)

  C1 = plt.contour(XX, YY, fr2d, levels = [0.0], linewidths = 2.0, colors = 'k')
  C2 = plt.contour(XX,YY, fi2d, levels = [0.0], linewidths = 2.0,colors = 'r')

# This is not needed but interesting to look at

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
              print('')
              print('S ', 'step ', 'm', 'k ', 'omega_r ', 'omega_i ', S, nk, m, k, ppx, ppy)
              print(px[maxindex],py[maxindex])

  nk = nk + 1
##^

#------------------------------------------------------------
# Write/plot the dispersion relationship
#------------------------------------------------------------

k_arr = np.linspace(0.001,kmax,nnk)

beta_arr = np.zeros(nnk)

# Write data to file:

with open('dispersion_relation_rankine_unstable.txt', 'w') as f:
  f.write('k (m-1)  omega_r (s-1)  omega_1 (s-1)\n')

  for i in np.arange(nnk):

    wstr = str(k_arr[i])+' '+str(or_arr[i])+' '+str(oi_arr[i])+'\n'
    f.write(wstr)
    print(k_arr[i], or_arr[i], oi_arr[i])
f.close()

plt.clf() # Somehow Python wants to plot the integrated previous figure

plt.plot(k_arr, or_arr, 'k', linewidth = 2.0)
plt.plot(k_arr, oi_arr, 'r', linewidth = 2.0)
#plt.savefig("./dispersion.png", dpi=120)
plt.show()

#------------------------------------------------------------
# The End.
#------------------------------------------------------------


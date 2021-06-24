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

# Solving the dispersion relation by Rotunno (1978, JFM).
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

m         = 2   # Azimuthal wavenumber 
Omega     = 1.0 # Base-state angular velocity
W         = 1.0 # Base-state axial velocity
R         = 1.0 # Cylinder radius
nmax      = 301 # Number of grid points in omega_r and omega_i directions
nnk       = 31  # Number of points on k axis of dispersion plot
kmax      = 5.0 # Maximum axial wavelength k
omega_max = 10.0 # Maximum wave frequency

# To make sure the solver only picks unstable modes, omega_i must be larger than eps

eps = 0.001 # Use a smaller number if nmax is increased

#-------------------------------------------------
# Start
#-------------------------------------------------

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

  omega_r = np.linspace(0, omega_max, nmax)
  omega_i = np.linspace(0, omega_max, nmax) * 1j

  X, Y = np.meshgrid(omega_r, omega_i)

  omega = X + Y

# Bessel function (1st kind) and its derivative

  I_m     = scipy.special.iv (m, k)
  I_mp    = scipy.special.ivp(m, k)

# Modified Bessel function (2nd kind) and its derivative

  K_m     = scipy.special.kv (m, k)
  K_mp    = scipy.special.kvp(m, k)
  
# Dispersion relation

  g1 = omega + W*k + 1.E-12   # inner (V = Omega = 0; W < 0)
  g2 = omega - m*Omega - W*k  # outer (W > 0)

  lhs = Omega**2 * R/(g1*g2) * I_mp * K_mp 
  rhs =-g2/(g1*k) * I_mp * K_m + g1/(g2*k) * K_mp * I_m

  D = lhs - rhs

  fr2d = D.real
  fi2d = D.imag

# Graphical solution (do not delete this part; C1 and C2 are used below)

  ox = np.linspace(0, omega_max, nmax)
  oy = np.linspace(0, omega_max, nmax)

  XX, YY = np.meshgrid(ox,oy)

  C1 = plt.contour(ox, oy, fr2d, levels = [0.0], linewidths = 2.0, colors = 'k')
  C2 = plt.contour(ox, oy, fi2d, levels = [0.0], linewidths = 2.0, colors = 'r')

# Not needed for the solved but interesting to look at

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
# Write/plot the dispersion relationship
#------------------------------------------------------------

k_arr = np.linspace(0,kmax,nnk)

# Write data to file:

with open('dispersion_relation_vortex_sheet.txt', 'w') as f:
  f.write('k (m-1)  omega_r (s-1)  omega_1 (s-1)\n')

  for i in np.arange(nnk):

    wstr = str(k_arr[i])+' '+str(or_arr[i])+' '+str(oi_arr[i])+'\n'
    f.write(wstr)
    print(k_arr[i], or_arr[i], oi_arr[i])
f.close()

# Plot preview

plt.clf() # Python wants to plot the integrated previous figure

plt.plot(k_arr, oi_arr, 'r', linewidth = 2.0)

plt.xlim(0, kmax)
plt.ylim(0,7.0)
plt.savefig("./dispersion.png", dpi=120)
plt.show()

#------------------------------------------------------------
# The End
#------------------------------------------------------------

# vortex_waves

This repository accompanies the review paper:

Dahl, J. M. L., 2021: Centrifugal Waves in Tornado-like Vortices: Kelvin's Solutions and Their Applications to Multiple-Vortex
Development and Vortex Breakdown. Mon. Wea. Rev., 149, 3173-3216, https://doi.org/10.1175/MWR-D-20-0426.1.

The repository contains scripts and a few text files used to generate the figures in the paper, as well as a document with clarifications and a list of errata. Note that many figures have been touched up by adding/modifying figure and axis labels (i.e., not all plots are precisely identical to the published version).

The following files are included:

1. clarifications.pdf: Some clarifications on the paper and corrections of minor typos.

2. couette_dispersion_m0.py: Python script that calculates and plots the dispersion relation for axisymmetric waves in a
   vortex bounded by cylindrical walls. Fig. 7 in the manuscript.
  
3. couette_dispersion_m1.py: Python script that calculates and plots the dispersion relation for spiral (m = 1) waves in a
   vortex bounded by cylindrical walls. Fig. 10 in the manuscript. 
   
4. dispersion_loiseleux_mm1_S1.txt: Text output of rankine_unstable_m1.py, used to plot Fig. 24 for m = -1

5. dispersion_loiseleux_mp1_S1.txt: Text output of rankine_unstable_m1.py, used to plot Fig. 24 for m = +1

6. dispersion_rotunno_m2_S1.txt: Text output of rotunno78.py, used to plot Fig. 26 for m = 2

7. dispersion_rotunno_m3_S1.txt: Text output of rotunno78.py, used to plot Fig. 26 for m = 3

8. hollow_vortex.py: Plots the dispersion relation for a hollow vortex; this scenario is not considered in the MWR paper (but it was considered by Kelvin).

9. rankine_axial_dispersion_stable.py: Python script that calculates and plots the dispersion relation for stable spiral (m = 1) waves in a
   Rankine vortex with an axial core jet.  This script also reads dispersion_loiseleux_mm1_S1.txt to plot the frequency
   of the unstable mode; Fig. 23 in the manuscript. 
    
10. rankine_pure_dispersion_cograde.py: Python script that calculates and plots the dispersion relation for axisymmetric (m = 0) and
   cograde spiral (m = 1) waves in a Rankine vortex (no axial flow). Figs. 16 and 18 in the manuscript. 
  
11. rankine_pure_dispersion_retrograde.py: Python script that calculates and plots the dispersion relation for retro-/countergrade spiral (m = 1) waves in a
   Rankine vortex (no axial flow). Fig. 19 in the manuscript. 
   
12. rankine_unstable_m1.py: Python script that calculates and writes out the dispersion relation of unstable spiral waves in the Rankine vortex
    with axial jet. Fig. 24 in the manuscript. (output: dispersion_loiseleux_mm1_S1.txt, dispersion_loiseleux_mp1_S1.txt)
  
13. rotunno78.py: Python script that calculates and writes out the dispersion relation of unstable spiral waves in the cylindrical vortex
    sheet with rising motion in its periphery and sinking motion in its interior. Fig. 26 in the manuscript. 
    (output: dispersion_rotunno_m2_S1.txt, dispersion_rotunno_m3_S1.txt)

14. spiral_waves_structure.py: Python script plotting the horizontal structure of the spiral modes in bounded and Rankine vortices; 
    Figs. 11, 12, 20, 21, 25, 27 in the manuscript.

15. structure_3d_loops.py: Python script that creates a large number of gifs that show 3D renderings of the spiral modes in the stable and
    unstable scenarios (supplemental material of the paper).  To create a loop, use e.g.: convert -loop 100 -delay 0.15 *png animation_NAME_OF_YOUR_CHOICE.gif

16. structure_meridional_m0.py: Python script that plots the meridional structure [(r,z) plane] of the axisymmetric modes in the bounded and Rankine
    vortex cases; Figs. 8, 9, 17 in the manuscript.

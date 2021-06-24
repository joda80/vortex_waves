# vortex_waves

This repository containes scripts and a few text files used to generate the figures in the MWR review article. Note that all
figure labels have been added/modified with Adobe Illustrator (i.e., not all plots are identical to the published version).

Dahl, J. M. L., 2021: Centrifugal Waves in Tornado-like Vortices: Kelvin's Solutions and Their Applications to Multiple-Vortex
Development and Vortex Breakdown. Mon. Wea. Rev., submitted.

The following files are included:

1. couette_dispersion_m0.py: Python script that calculates and plots the dispersion relation for axisymmetric waves in a
   vortex bounded by cylindrical walls. Fig. 7 in the manuscript.
  
2. couette_dispersion_m1.py: Python script that calculates and plots the dispersion relation for spiral (m = 1) waves in a
   vortex bounded by cylindrical walls. Fig. 10 in the manuscript. 
   
3. dispersion_loiseleux_mm1_S1.txt: Text output of rankine_unstable_m1.py, used to plot Fig. 24 for m = -1

4. dispersion_loiseleux_mp1_S1.txt: Text output of rankine_unstable_m1.py, used to plot Fig. 24 for m = +1

5. dispersion_rotunno_m2_S1.txt: Text output of rotunno78.py, used to plot Fig. 26 for m = 2

6. dispersion_rotunno_m3_S1.txt: Text output of rotunno78.py, used to plot Fig. 26 for m = 3

7. rankine_axial_dispersion_stable.py: Python script that calculates and plots the dispersion relation for stable spiral (m = 1) waves in a
   Rankine vortex with an axial core jet.  This script also reads dispersion_loiseleux_mm1_S1.txt to plot the frequency
   of the unstable mode; Fig. 23 in the manuscript. 
    
8. rankine_pure_dispersion_cograde.py: Python script that calculates and plots the dispersion relation for axisymmetric (m = 0) and
   cograde spiral (m = 1) waves in a Rankine vortex (no axial flow). Figs. 16 and 18 in the manuscript. 
  
9. rankine_pure_dispersion_retrograde.py: Python script that calculates and plots the dispersion relation for retro-/countergrade spiral (m = 1) waves in a
   Rankine vortex (no axial flow). Fig. 19 in the manuscript. 
   
10. rankine_unstable_m1.py: Python script that calculates and writes out the dispersion relation of unstable spiral waves in the Rankine vortex
    with axial jet. Fig. 24 in the manuscript. (output: dispersion_loiseleux_mm1_S1.txt, dispersion_loiseleux_mp1_S1.txt)
  
11. rotunno78.py: Python script that calculates and writes out the dispersion relation of unstable spiral waves in the cylindrical vortex
    sheet with rising motion in its periphery and sinking motion in its interior. Fig. 26 in the manuscript. 
    (output: dispersion_rotunno_m2_S1.txt, dispersion_rotunno_m3_S1.txt)

12. spiral_waves_structure.py: Python script plotting the horizontal structure of the spiral modes in bounded and Rankine vortices; 
    Figs. 11, 12, 20, 21, 25, 27 in the manuscript.

13. structure_3d_loops.py: Python script that creates a large number of gifs that show 3D renderings of the spiral modes in the stable and
    unstable scenarios (supplemental material of the paper).  To create a loop, use e.g.: convert -loop 100 -delay 0.15 *png animation_NAME_OF_YOUR_CHOICE.gif

14. structure_meridional_m0.py: Python script that plots the meridional structure [(r,z) plane] of the axisymmetric modes in the bounded and Rankine
    vortex cases; Figs. 8, 9, 17 in the manuscript.

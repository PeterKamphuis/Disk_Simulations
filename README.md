# Disk_Simulations
A fortran program to create artificial fits images of a exponentional disk.

This code is used in Kamphuis et al. 2007 (https://ui.adsabs.harvard.edu/abs/2007A%26A...468..951K/abstract) and Voigtlander et al 2013 (https://ui.adsabs.harvard.edu/abs/2013A%26A...554A.133V/abstract).



The fortran version compiles with
gfortran -o Disk_Simulation Exponential_Disk_Simulation.f -L/usr/lib/ -lcfitsio

Or replace /usr/lib/ with whereever your cfitsio is located. The code can then be run with:

./Disk_Simulation Input_Parameters.txt

Where Input_Parameters is the name of the file containing the input parameters (see below).

It compiles without errors and runs but is rather slow
Seems to produce the correct output though.


The parameter list clearly needs some cleanup though.
Here is the explanation of the Parameters. These are read in order including the comment lines so the input has to fixed.

# File name for output disk
File name for output disk
#  Inclination of disk
Input Inclination for the model
#  Disk scale length (kpc)
Scale length of the Disk
#  Disk scale height (kpc)
Scale height of the disk
#  Dust scale length (kpc)
Scale length of the dust disk
#  Dust scale height (kpc)
Scale height of the dusk disk
#  Optical thickness of disk
Opacity of the dust disk
#  Bulge to disk ratio
This the factor the intensity in the bulge is multiplied with before being combined with the disk.
#  B0
Scale length of the bulge, the bulge is treated in exactly the same way as the disk, i.e. an exponentional but with the bulge parameters
#  Ellipticity of bulge
Factor to multiple the bulge scale length with to get bulge scale height
#  Logarithmic binning of disk
If True the lower x and z boundaries are set to 1 if < 0 and else to log10(boun) the upper boundaries are set to log10(10)
Presumably to get log increasing sizes but this is not well tested.
#  Fold disk to make disk symmetric
Option to caluclate only a quarter of the disk and reuse this value for all four quarters of the disk.
# Truncation Radius (kpc)
Truncation radius of the galaxy, i.e. all disk. Beyond this radius in the galaxy all emission and absorption is set to 0.
#  X lower limit
lower x limit of the fits image in kpc
#  X upper limit
upper x limit of  the fits image in kpc
#  Z lower limit
lower y limit of  the fits image in kpc
#  Z upper limit
upper y limit of  the fits image in kpc
#  Number of pixels in x
size in x, in pixels
#  Number of pixels in z
size in y, in pixels
#  Number of steps in Gauss-legendre integration
Number of steps in Gauss-legendre integration, the higher number the higher the accuracy of the model.
# Inner Hole radius (kpc)
between these radii all emission is set to 0
# Outer Hole radius (kpc)
between these radii all emission is set to 0
#  Boundary 1
The inner boundary of the dust disk, at radii smaller no absorption is calculated.
#  Boundary 2
The outer boundary of the dust disk, at larger radii no absorption is calculated.
#  Boundary 3
The inner boundary of a secondary dust disk
#  Boundary 4
The inner boundary of a secondary dust disk
# Distance(Mpc)
Used for the correct pixel scale in arcsec in fits file.
# Create a cube
If yes a 3D velocity cube based on the rc is created else only an image of the disk is created and the following parameters are not required.
#  File name for output Cube
File name for output Cube fits file
#  File with rotation curve
A file name should always be present
# lag in km/s/kpc
decline of RC as function of height above the disk
#  Velocity dispersion (isotropic)(central)
Inner velocity dispersion
#  Velocity dispersion (isotropic)(outer)
Outer velocity dispersion, the code interpolates linearly from inner to outer dispersion.
#  Velocity dispersion of bulge
Velocity dispersion of the bulge
#  Starting velocity (km/s)
Minimum velocity away from the central velocity
#  Step in velocity (km/s)
Channel width
#  Central velocity of galaxy (km/s)
systemic velocity
#  Number of pixels in v
size in z, in pixels

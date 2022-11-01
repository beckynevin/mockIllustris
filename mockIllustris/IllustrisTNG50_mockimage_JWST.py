# Nov 1 2022 - updating to be a .py file
### Goal: Create a mock image from a TNG50 SKIRT datacube
#  This notebook will walk through the steps to go from the SKIRT output cube into a reasonable mock of the user's choice of filter (from JWST or HST). 

#  The user can also set the redshift of the image and the size of the image cutout.

#  Things this notebook does:
#  1. Imports SED from SKIRT and describes the file (from Xuejian Shen).
#  2. Applies cosmological dimming and IGM absorption (from Xuejian Shen).
#  3. Uses SEDPY to create an image of a given filter.
#  4. Applies an appropriate PSF.
#  5. Rebins to the pixelscale of the instrument.
#  6. Introduces appropriate background residual.

#  Things I haven't done yet: 
#  1. Centering software to center on the galaxy of choice.
#  2. Determine (if we need) and implement background galaxies.
#  3. Potentially create an error image.

#  Things that need to be double checked:
#  1. The cosmological dimming.
#  2. Adding the background residuals in the appropriate units.

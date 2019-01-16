# takes saved plot images and creates gif
#NOTE: on scinet, need to use python 2.7.8 for this script to work -- tobytes function doesn't exist in python 2.7.5 (numpy 1.7)
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import colormaps as cmaps
import os.path
import imageio

zInit=0.0
zFinal=4.6
n=22
dz=0.2
i=1
nside = 2048

gifArray = []

while i<n:
    print "i:", i
    zMin = i*dz
    zMax = (i+1)*dz
    zMinCib = (i+1)*dz
    zMaxCib = (i+2)*dz
    
    imgFilename = "/scratch2/r/rbond/phamloui/pics/jul17_lensed_fwhm_0.0035/lensed_zmin%.1f_zmax%.1f.png" % (zMinCib, zMaxCib)
    imgSlice = imageio.imread(imgFilename)
    print imgSlice.shape
    gifArray.append(imgSlice)

    i += 1

imageio.mimsave("/scratch2/r/rbond/phamloui/pics/jul17_lensed_fwhm_0.0035/lensed_anim.gif", gifArray, duration=0.7, subrectangles=True)

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from astropy.io import fits
import os
from scipy import interpolate

zInit=0.0
zFinal=4.6
n=22
dz=0.2
i=0

z = np.arange(zInit+(2*dz), zFinal+dz, dz)
nside = 2048
lmax = 3 * nside - 1

f1 = plt.figure()
plt.title("CIB Power Spectra @ l=500 per shell")
plt.xlabel(r'$z$', fontsize=40)
plt.ylabel(r"$dC_l/dz$", fontsize=40)
print z
print z.shape[0]
cl = []
while i<n:
    print "i:", i
    # if i%3 != 0:
    # i += 1
    # continue

    zMin = i*dz
    zMax = (i+1)*dz
    zMinCib = (i+1)*dz
    zMaxCib = (i+2)*dz
    lensedClFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_lensed_cl/lensed_cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_ns2048_tot.dat" % (zMinCib, zMaxCib)
    lensedCl = hp.read_cl(lensedClFilename)
    pt = lensedCl[500]
    # plt.plot(z[i], pt, 'r-')
    cl.append(pt)
    i += 1

z_interp = np.arange(0,4.7, 0.1)
cl_interp_func = interpolate.splrep(z, cl, s=0)
cl_interp = interpolate.splev(z_interp, cl_interp_func, der=0) 

plt.plot(z, cl, 'b-', linewidth=2)
plt.show()

# plots gradphi over phi 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
nside=2048
pixels = np.arange(0,nside**2*12)
radToDeg = 180./np.pi
theta, phi = hp.pix2ang(nside, pixels)
mollProj = hp.projector.MollweideProj()
x,y = mollProj.ang2xy(theta, phi)
theta = 90 - np.degrees(theta) #lat
phi = np.degrees(phi) #lon

seismic_cmap = cm.get_cmap('seismic')
seismic_cmap.set_under('w')

gradphiR, gradphiI = hp.read_map("/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_phi/8Gpc_n2048_nb18_nt16_phi_sis_2_ns2048_zmin0.0_zmax2.8_hp_grad.fits", field=(0,1))
# gradphiR = hp.ud_grade(gradphiR, 512)
# gradphiI = hp.ud_grade(gradphiI, 512)
gravPotential = hp.alm2map(hp.read_alm("/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_phi/8Gpc_n2048_nb18_nt16_phi_sis_2_ns2048_zmin0.0_zmax2.8_hp.fits"), 2048)
unlensed = hp.read_map("/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_unlensed/cib_fullsky_ns2048_zmin2.80_zmax3.00_nu217_ns2048_tot.fits")
lensed = hp.read_map("/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_lensed/lensed_cib_fullsky_ns2048_zmin2.80_zmax3.00_nu217_ns2048_tot.fits")
# print lensed[np.isnan(lensed)].shape
# gradIndices = np.random.choice(nside**2*12, nside**2*12/40000, replace=False)
# gradIndices = pixels[((phi<=90) | (phi>=270)) & (np.absolute(theta)<=45)]
# gradIndices = np.random.choice(gradIndices, gradIndices.shape[0]/400, replace=False)

# get nan coordinates and overplot there
lensedUnseen = np.copy(lensed)
lensedNanToZero = np.copy(lensed)
nanPos = np.isnan(lensed)
nanIndices = np.where(nanPos)
# print nanIndices
lensedNanToZero[nanPos] = 0
lensedUnseen[nanPos] = hp.UNSEEN
cibMin = np.minimum(unlensed.min(), lensedNanToZero.min())
cibMax = np.maximum(unlensed.max(), lensedNanToZero.max())

fig=plt.figure(1)
hp.mollview(gravPotential, fig=1, xsize=2048, hold=True)
plt.quiver(x[nanIndices], y[nanIndices], gradphiI[nanIndices], gradphiR[nanIndices], units='x', width=0.0006, scale=0.045, pivot='mid') #lower scale -> larger arrows

plt.figure(2)
hp.mollzoom(unlensed, fig=2, xsize=2048, title="Unlensed", min=cibMin, max=cibMax)
hp.set_g_clim(cibMin, cibMax)
plt.figure(3)
hp.mollzoom(lensedUnseen, fig=3, xsize=2048, title="Lensed", min=cibMin, max=cibMax)
hp.set_g_clim(cibMin, cibMax)

# plt.figure(4)
# diffMap = lensed - unlensed
# diffRange = np.absolute(diffMap).max()
# diffScheme = colors.SymLogNorm(linthresh=1.2e-20, linscale=0.03, vmin=-diffRange, vmax=diffRange) # bigger linthresh -> less stuff
# diffProj = hp.mollview(diffMap, fig=4, xsize=2048, hold=True, return_projected_map=True, cmap=seismic_cmap, min=-diffRange, max=diffRange, norm=diffScheme)
# plt.quiver(x[gradIndices], y[gradIndices], gradphiI[gradIndices], gradphiR[gradIndices], units='width', width=0.0006, scale=0.10, pivot='mid')

plt.show()

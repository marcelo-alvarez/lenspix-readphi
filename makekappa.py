# combine marcelo's z<4.5 map with GRF from CAMB for remaining 4.5<z<1100
import numpy as np
import healpy as hp
from scipy import interpolate
import matplotlib.pyplot as plt

lmax=4000
nside=2048
CMB_outputscale = 7.4311e12 #straight outta CAMB (default value)

# load CAMB z<1100
#camb_ell, camb_C_phi = np.loadtxt("/scratch2/r/rbond/engelen/lenspower/CAMB/peakpatch_lenspotenticalCls.dat", usecols=(0,5), unpack=True)
camb_ell, camb_C_phi = np.loadtxt("/scratch2/r/rbond/engelen/lenspower/CAMB/peakpatch_nonlinear0_scalCls.dat", usecols=(0,4), unpack=True)
camb_ell = np.insert(camb_ell, (0,0), (0,1)) #missing ell=0,1, add zeros
camb_C_phi = np.insert(camb_C_phi, (0,0), (0,0))

#camb_CL_phi = camb_C_phi * 2*np.pi / (camb_ell*(camb_ell+1))**2 #camb's lenspotentialCls.dat actually stores C_dd = [l(l+1)]^2 * C_L^phi / 2pi
camb_CL_phi = camb_C_phi/ (camb_ell**4) / CMB_outputscale #scalCls factor
camb_CL_kappa = camb_CL_phi * (camb_ell*(camb_ell+1)/2.0)**2
camb_CL_kappa[0] = 0
camb_CL_kappa[1] = 0
#print "linear camb kappa:", camb_CL_kappa[0:5]

# load CAMB z<1100 + halofit
camb_ell_nonlin, camb_C_phi_nonlin = np.loadtxt("/scratch2/r/rbond/engelen/lenspower/CAMB/peakpatch_nonlinear3_scalCls.dat", usecols=(0,4), unpack=True)
camb_ell_nonlin = np.insert(camb_ell_nonlin, (0,0), (0,1)) #missing ell=0,1, add zeros
camb_C_phi_nonlin = np.insert(camb_C_phi_nonlin, (0,0), (0,0))

camb_CL_phi_nonlin = camb_C_phi_nonlin / (camb_ell_nonlin**4) / CMB_outputscale
camb_CL_kappa_nonlin = camb_CL_phi_nonlin * (camb_ell_nonlin*(camb_ell_nonlin+1)/2.0)**2
camb_CL_kappa_nonlin[0] = 0
camb_CL_kappa_nonlin[1] = 0
#print "nonlinear camb kappa:", camb_CL_kappa_nonlin[0:5]

# output the 'camb w/ halofit all-z' as kappa cl and map
#camb_CL_kappa_nonlin_map = hp.synfast(camb_CL_kappa_nonlin, nside, lmax=lmax)
#print "writing full z kappa..."
#hp.write_map("/scratch2/r/rbond/phamloui/lenspix_files/kappa_maps/all_z_kappa_with_halofit.fits", camb_CL_kappa_nonlin_map)
#hp.write_cl("/scratch2/r/rbond/phamloui/lenspix_files/all_z_kappa_with_halofit.dat", camb_CL_kappa_nonlin[0:lmax+1])

# load camb z<4.5
alex_ell_raw, alex_CL_kappa_raw = np.loadtxt("/scratch2/r/rbond/engelen/peakpatch/data/cl_kappa_restricted_nonlinear0.txt", usecols=(0,1), unpack=True)
alex_interp_func = interpolate.splrep(alex_ell_raw, alex_CL_kappa_raw, s=0)
alex_ell_interp = np.arange(0,lmax+1)
alex_CL_kappa_interp = interpolate.splev(alex_ell_interp, alex_interp_func, der=0)
alex_CL_kappa_interp[0] = 0
alex_CL_kappa_interp[1] = 0
alex_CL_phi_interp = alex_CL_kappa_interp / (alex_ell_interp * (alex_ell_interp+1) / 2.0)**2
#print "linear alex kappa:", alex_CL_kappa_interp[0:5]

# load camb z<4.5 + halofit
alex_ell_raw_nonlin, alex_CL_kappa_raw_nonlin = np.loadtxt("/scratch2/r/rbond/engelen/peakpatch/data/cl_kappa_restricted_nonlinear3.txt", usecols=(0,1), unpack=True)
alex_interp_func_nonlin = interpolate.splrep(alex_ell_raw_nonlin, alex_CL_kappa_raw_nonlin, s=0)
alex_ell_interp_nonlin = np.arange(0,lmax+1)
alex_CL_kappa_interp_nonlin = interpolate.splev(alex_ell_interp_nonlin, alex_interp_func_nonlin, der=0)
alex_CL_kappa_interp_nonlin[0] = 0
alex_CL_kappa_interp_nonlin[1] = 0
alex_CL_phi_interp_nonlin = alex_CL_kappa_interp_nonlin / (alex_ell_interp_nonlin * (alex_ell_interp_nonlin+1) / 2.0)**2
#print "nonlinear alex kappa: ", alex_CL_kappa_interp_nonlin[0:5]

# load peakpatch z<4.5
marcelo_kappa_map = hp.read_map("/scratch2/r/rbond/malvarez/peakpatch/lightcones/octant/8Gpc_n4096_nb18_nt16/kappa/8Gpc_n4096_nb18_nt16_kap_fph.fits")
marcelo_CL_kappa = hp.anafast(marcelo_kappa_map, lmax=lmax)

remainder_CL_kappa = camb_CL_kappa[0:lmax+1] - alex_CL_kappa_interp
remainder_CL_phi = camb_CL_phi[0:lmax+1] - alex_CL_phi_interp
#print "linear remainder kappa:", remainder_CL_kappa[0:5]
remainder_kappa_map = hp.synfast(remainder_CL_kappa, nside)
remainder_CL_kappa_grf = hp.anafast(remainder_kappa_map, lmax=lmax) #just remainder_CL_kappa with GRF

remainder_CL_kappa_nonlin = camb_CL_kappa_nonlin[0:lmax+1] - alex_CL_kappa_interp_nonlin
remainder_CL_phi_nonlin = camb_CL_phi_nonlin[0:lmax+1] - alex_CL_phi_interp_nonlin
#print "nonlinear remainder kappa:", remainder_CL_kappa[0:5]
remainder_kappa_map_nonlin = hp.synfast(remainder_CL_kappa_nonlin, nside)
remainder_CL_kappa_grf_nonlin = hp.anafast(remainder_kappa_map_nonlin, lmax=lmax)

total_CL_kappa = marcelo_CL_kappa + remainder_CL_kappa_grf # want to set ell=0,1 to zero before getting total kappa map
total_CL_kappa[0]=0
total_CL_kappa[1]=0
total_kappa_map = hp.synfast(total_CL_kappa, nside, lmax=lmax)

total_CL_kappa_nonlin = marcelo_CL_kappa + remainder_CL_kappa_grf_nonlin
total_CL_kappa_nonlin[0] = 0
total_CL_kappa_nonlin[1] = 0
total_kappa_map_nonlin = hp.synfast(total_CL_kappa_nonlin, nside, lmax=lmax)

#plt.figure()
#plt.xlabel(r"$l$", fontsize=30)
#plt.ylabel(r"$C_L^{\phi\phi}$", fontsize=30)
#plt.loglog(camb_ell[0:lmax+1], camb_CL_phi[0:lmax+1], 'b--', label="camb 0->1100")
#plt.loglog(alex_ell_interp, alex_CL_phi_interp, 'r--', label="alex 0->4.5")
#plt.loglog(alex_ell_interp, remainder_CL_phi, 'g--', label="phi remainder 4.5->1100")
#legend3 = plt.legend(loc="lower left", shadow=True)
#frame3 = legend3.get_frame()
#frame3.set_facecolor('0.90')

#l(l+1)/2pi
#camb_CL_kappa *= camb_ell*(camb_ell+1)/2/np.pi
#alex_CL_kappa_interp *= alex_ell_interp*(alex_ell_interp+1)/2/np.pi
#remainder_CL_kappa *= alex_ell_interp*(alex_ell_interp+1)/2/np.pi

#camb_CL_kappa_nonlin *= camb_ell*(camb_ell+1)/2/np.pi
#alex_CL_kappa_interp_nonlin *= alex_ell_interp_nonlin*(alex_ell_interp_nonlin+1)/2/np.pi
#remainder_CL_kappa_nonlin *= alex_ell_interp_nonlin*(alex_ell_interp_nonlin+1)/2/np.pi

#total_CL_kappa *= alex_ell_interp*(alex_ell_interp+1)/2/np.pi
#total_CL_kappa_nonlin *= alex_ell_interp*(alex_ell_interp+1)/2/np.pi
#marcelo_CL_kappa *= alex_ell_interp*(alex_ell_interp+1)/2/np.pi

plt.figure()
plt.xlabel(r"$l$", fontsize=30)
plt.ylabel(r"$C_l^{\kappa\kappa}$", fontsize=30)
#plt.plot(alex_ell_raw, alex_CL_kappa_raw,'rx', label="alex 0->4.5")
plt.loglog(camb_ell[0:lmax+1], camb_CL_kappa[0:lmax+1], 'b--', label="camb 0->1100") 
plt.loglog(alex_ell_interp, alex_CL_kappa_interp, 'r--', label="alex 0->4.5 interpolated") #all spectra ell should be the same
plt.loglog(alex_ell_interp, remainder_CL_kappa, 'g--', label="remainder 4.5->1100")

plt.loglog(camb_ell[0:lmax+1], camb_CL_kappa_nonlin[0:lmax+1], 'm:', label="camb nonlinear 0->1100")
plt.loglog(alex_ell_interp_nonlin, alex_CL_kappa_interp_nonlin, 'y:', label="alex nonlinear 0->4.5")
plt.loglog(alex_ell_interp_nonlin, remainder_CL_kappa_nonlin, 'k:', label="nonlinear remainder 4.5->1100")
plt.loglog(alex_ell_interp, remainder_CL_kappa_grf, 'm', label="new map 4.5->1100")
plt.loglog(alex_ell_interp, total_CL_kappa, 'k', label="marcelo map+new map")
plt.loglog(alex_ell_interp, total_CL_kappa_nonlin, 'c', label="marcelo map+new nonlinear map")

plt.loglog(alex_ell_interp, marcelo_CL_kappa, '0.75', label="marcelo map")

legend2 = plt.legend(loc="lower right", shadow=True)
frame2 = legend2.get_frame()
frame2.set_facecolor('0.90')

plt.show() 

#print "writing kappa..."
#day_root = "jun1_"
#hp.write_map("/scratch2/r/rbond/phamloui/lenspix_files/kappa_maps/"+day_root+"nonlinear0_kappa_for_julian_cmb.fits", total_kappa_map)
#hp.write_map("/scratch2/r/rbond/phamloui/lenspix_files/kappa_maps/"+day_root+"nonlinear3_kappa_for_julian_cmb_only_halofit.fits", total_kappa_map_nonlin)

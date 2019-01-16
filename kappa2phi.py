import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

def kap2phi(field_kappa, halo_kappa, unlensed_primary, phi_alm_file, writeMap=False, phi_map_file=None, lens_lmax=None):
    print "===============kappa2phi==============="
    nside = hp.get_nside(field_kappa)    
    lmax = hp.Alm.getlmax(len(unlensed_primary))
    if not lens_lmax:
        lens_lmax = lmax
    print "LMAX:", lmax
    #print "----Done loading maps."

    #combine
    print "----Combining field and halo kappa..."
    kappa_map = field_kappa + halo_kappa
    print "----Done combining."
    
    #convert to alm
    print "----Converting kappa map to alm..."
    kappa_lm = hp.map2alm(kappa_map, lmax=lens_lmax)
    print "----Done."
    
    #convert to phi (grav potential)
    print "----Converting kappa to phi..."
    l,m = hp.Alm.getlm(lens_lmax)
    phi_lm = kappa_lm * (2.0 / (l*(l+1.0)))
    phi_lm[l==0] = 0

    print "----Writing phi alm to file..."
    hp.write_alm(phi_alm_file, phi_lm)
    print "----Done."

    if writeMap and phi_map_file:
        #convert to map
        print "----Converting phi to map..."
        phi_map = hp.alm2map(phi_lm, nside, lmax=lmax)
        print "----Writing phi map to file..."
        hp.write_map(phi_map_file, phi_map)
        print "----Done."

    print "=============kappa2phi end============="
    return lmax

if __name__ == "__main__":
    #todo: file selection
    #input filenames
    #field_kappa_file = "kappa_maps/8Gpc_n4096_nb23_nt18_kap_field.fits"
    #halo_kappa_file = "kappa_maps/8Gpc_n4096_nb23_nt18_kap_halo.fits"
    #unlensed_primary_file = "FromNERSC/ffp10_unlensed_scl_cmb_000_alm.fits"
    #output filenames
    #comb_kappa_file = "kappa_maps/8Gpc_n4096_nb23_nt18_kap_comb.fits"
    #phi_alm_file = "kappa_maps/comb_z_4.6_n4096_lmax5120_phi_alm.fits"
    #phi_map_file = "kappa_maps/comb_z_4.6_n4096_lmax5120_phi_map.fits"
    #kap2phi(field_kappa_file, halo_kappa_file, unlensed_primary_file, comb_kappa_file, phi_alm_file, writeMap=True, phi_map_file=phi_map_file)
    pass

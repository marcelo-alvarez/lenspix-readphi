#!/bin/bash
# MOAB/Torque submission script for SciNet GPC
#
#PBS -l nodes=5:ppn=8,walltime=3:00:00
#PBS -N lensing_maps_slicing

# DIRECTORY TO RUN - $PBS_O_WORKDIR is directory job was submitted from
cd $PBS_O_WORKDIR

zi=0.0
zf=4.6
nz=22
dz=0.2
i=1

while [ $i -lt $nz ]; do
    echo $i
    zmin=`awk "BEGIN {print $zi+$i*$dz}"`
    zmax=`awk "BEGIN {print $zi+($i+1)*$dz}"`

    zmincib=`awk "BEGIN {print $zi+($i+1)*$dz}"` #next redshift from kappa
    zmaxcib=`awk "BEGIN {print $zi+($i+2)*$dz}"` #ie kappa of 0<z<0.2 lenses cib 0.2<z<0.4

    #round zcib to fit with cib filenames
    zmincib=`printf "%0.2f\n" $zmincib`
    zmaxcib=`printf "%0.2f\n" $zmaxcib`

    #lens cib z slice with kappa map of structure integrated up to z slice - kappa converted to phi prior to lensing
    # <kappamap_in> <cibmap_in> <phimap_out> <cibmap_out>
    # python lens.py /scratch2/r/rbond/phamloui/lenspix_files/lensing/kappa_z4pt6_z0.0_z${zmax}_nside2048_hp.fits /scratch2/r/rbond/phamloui/lenspix_files/cib/cib_fullsky_ns2048_zmin${zmincib}_zmax${zmaxcib}_nu217_tot.fits /scratch2/r/rbond/phamloui/lenspix_files/lensing/phi_z4pt6_z0.0_z${zmax}_nside2048_hp.fits /scratch2/r/rbond/phamloui/lenspix_files/cib/lensed_cib_fullsky_ns2048_zmin${zmincib}_zmax${zmaxcib}_nu217_tot.fits -np 24

    python lens.py /scratch2/r/rbond/phamloui/lenspix_files/cib_kappa_zmin0.2/8Gpc_n4096_nb18_nt16_kap_sis_2_ns2048_zmin0.2_zmax${zmax}_hp.fits /scratch2/r/rbond/phamloui/lenspix_files/cib_v2_unlensed/cib_fullsky_ns2048_zmin${zmincib}_zmax${zmaxcib}_nu217_ns2048_tot_fwhm_0.0035.fits /scratch2/r/rbond/phamloui/lenspix_files/cib_phi_zmin0.2/8Gpc_n4096_nb18_nt16_phi_sis_2_ns2048_zmin0.2_zmax${zmax}_hp.fits /scratch2/r/rbond/phamloui/lenspix_files/cib_lensed_kap0.2/jul13_cib_fullsky_ns2048_zmin${zmincib}_zmax${zmaxcib}_nu217_ns2048_tot_fwhm_0.0035.fits -np 40
    i=$(($i+1))

done


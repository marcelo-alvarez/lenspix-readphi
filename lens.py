from kappa2phi import kap2phi
import numpy as np
import healpy as hp
from astropy.io import fits
import subprocess
import argparse
import os.path

def checkFitsType(fileName):
    '''
    Check if FITS file is a map or alm.
    '''
    hdulist = fits.open(fileName)
    try:
        nside = hdulist[1].header['NSIDE']
        fitsType = 'map'
    except KeyError:
        #assume alm
        fitsType = 'alm'
    except Exception:
        raise
    hdulist.close()
    return fitsType

def appendToFilename(filename, newStr):
    return "{0}_{2}.{1}".format(*filename.rsplit('.', 1) + [newStr])

#at simlens step, will need a phi and primary alm -- check inputs here
parser = argparse.ArgumentParser(description='Lens CMB primary.')
parser.add_argument('input_field_kappa', help='filename of the field kappa file (pref. HEALPix map)')
#parser.add_argument('input_halo_kappa', help='filename of the halo kappa file (HEALPix map)')
parser.add_argument('input_primary', help='filename of the unlensed primary file (pref. HEALPix alm)')
parser.add_argument('output_phi', help='filename of the outputted phi file (HEALPix alm)')
parser.add_argument('output_lensed', help='filename of the outputted lensed file (HEALPix map)')
parser.add_argument('-ol', '--output_lmax', help='new lmax for the lensed map')
parser.add_argument('-np', '--num_processes', help='-np from job submission')
args = parser.parse_args()

genericParams = 'bin/generic_params.ini'
specificParams = 'bin/specific_params.ini' #NOTE: cant save to home directory when running this script as a job
outFileRoot = "lensing_output_root" #used mainly for power spectra files - in the future just append _power to filenames from arguments?

if args.num_processes:
    numProc = args.num_processes
else:
    numProc = '1'

#input filenames
fieldKappaFile = args.input_field_kappa
#haloKappaFile = args.input_halo_kappa
unlensedPrimaryFile = args.input_primary

#output filenames                                                                                    
phiAlmFile = args.output_phi
lensedFile = args.output_lensed
gradPhiFile = appendToFilename(phiAlmFile, 'grad')

#get nside
hdulist = fits.open(fieldKappaFile)
nside = hdulist[1].header['NSIDE']
hdulist.close()

#no halo map yet, use zero map for now
zeroMapFile = "test_maps/zeros_%s.fits" % (nside)
haloKappaFile = zeroMapFile

if not os.path.exists(zeroMapFile):
    print "Creating zero map..."
    zeroMap = np.zeros(12 * (nside**2)) #create a zero map and use that in place of halo kappa, since we only have field kappa now
    hp.write_map(zeroMapFile, zeroMap)

#load maps
print "Loading maps..."
fieldKappaType = checkFitsType(fieldKappaFile)
primaryType = checkFitsType(unlensedPrimaryFile)
print "Loading field kappa..."
fieldKappa = hp.read_map(fieldKappaFile)
print "Loading halo kappa..."
haloKappa = hp.read_map(haloKappaFile)
if primaryType != 'alm':
    print "Primary is map - converting to alm..."
    unlensedPrimaryMap = hp.read_map(unlensedPrimaryFile)
    print "get filename"
    unlensedPrimaryFile = appendToFilename(unlensedPrimaryFile, "alm")
    print unlensedPrimaryFile
    print "convert"
    unlensedPrimary = hp.map2alm(unlensedPrimaryMap)
    print "write"
    hp.write_alm(unlensedPrimaryFile, unlensedPrimary)
    print "Saved new primary alm to", unlensedPrimaryFile
else:
    print "Loading primary..."
    unlensedPrimary = hp.read_alm(unlensedPrimaryFile)
unlenLmax = hp.Alm.getlmax(unlensedPrimary.shape[0])

#get lmax and create phi alm
if not os.path.exists(phiAlmFile):
    print "phi doesn't exist"
    phiLmax = kap2phi(fieldKappa, haloKappa, unlensedPrimary, phiAlmFile)
    #lmax = kap2phi(fieldKappaFile, haloKappaFile, unlensedPrimaryFile, phiAlmFile)
else:
    print "phi exists"
    phiAlm = hp.read_alm(phiAlmFile)
    phiLmax = hp.Alm.getlmax(phiAlm.shape[0])

lmax = max(unlenLmax, phiLmax)

print 'Obtained parameters NSIDE:', nside, 'and LMAX:', lmax

print 'Creating params file...'
if args.output_lmax:
    outputLmax = args.output_lmax
    print 'Output lmax set to', str(outputLmax) + '.'
else:
    outputLmax = lmax
    print 'Output lmax defaulted to input lmax', str(outputLmax) + '.'

subprocess.call(['cp', genericParams, specificParams])
#commas are used as delimiters for sed commands to avoid escaping (back)slashes -- this assumes there's no commas in the filenames
subprocess.call(['sed', '-i', 's,__NSIDEREPLACE__,' + str(nside) + ',g', specificParams])
subprocess.call(['sed', '-i', 's,__LMAXREPLACE__,' + str(lmax) + ',g', specificParams])
subprocess.call(['sed', '-i', 's,__OUTFILEROOTREPLACE__,' + outFileRoot + ',g', specificParams])
subprocess.call(['sed', '-i', 's,__OUTPUTLMAXREPLACE__,' + str(outputLmax) + ',g', specificParams])
subprocess.call(['sed', '-i', 's,__PRIMARYFILEREPLACE__,' + unlensedPrimaryFile + ',g', specificParams])
subprocess.call(['sed', '-i', 's,__PHIFILEREPLACE__,' + phiAlmFile + ',g', specificParams])
subprocess.call(['sed', '-i', 's,__GRADPHIFILEREPLACE__,' + gradPhiFile + ',g', specificParams])
subprocess.call(['sed', '-i', 's,__LENSEDFILEREPLACE__,' + lensedFile + ',g', specificParams])
print 'Params file created.'

print 'Running simlens...'
print 'Running "mpirun -np %s ./bin/simlens %s"' % (numProc, specificParams)
subprocess.call(['mpirun', '-np', str(numProc), './bin/simlens', specificParams])
print 'Simlens complete.'

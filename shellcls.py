import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from astropy.io import fits
import os
from pycamb_scripts.savitzky_golay import savitzky_golay

zInit=0.0
zFinal=4.6
n=22
dz=0.2
i=0

dPlot = 1

nside = 2048
lmax = 3 * nside - 1
summedCl = np.zeros(lmax+1)
summedLensedCl = np.copy(summedCl)
# colours = ['#ff7f00','#fdbf6f','#e31a1c','#fb9a99','#33a02c','#b2df8a','#1f78b4','#a6cee3']  
# colours = ['#e52d2d', '#e9542c', '#ed712a', '#ef8928', '#c0ce00', '#a7c304', '#8eb708', '#74ac0b', '#008910', '#1d9747', '#22a471', '#19b094', '#00e5d2', '#3bb2c5', '#4188b9', '#3a58ab', '#000fe5', '#7f02d7', '#b500c9', '#e00cbb', '#8e0047', '#480f29', '#190514'] #rainbow colour scheme
# colours = ['#1eff4f', '#23ec4a', '#26d946', '#28c942', '#29b93e', '#29a93a', '#289d36', '#289033', '#278430', '#267b2d', '#25712a', '#246828', '#225f25', '#215623', '#1f4d20', '#1d441d', '#1b3c1b', '#193418', '#162b15', '#142412', '#111c0f', '#0d170b', '#081007'] #green sequential
# colours = ['#ffbd8e', '#efb185', '#e1a77e', '#d49e77', '#c79470', '#ba8a69', '#ad8162', '#a37a5d', '#997358', '#8f6b52', '#86644d', '#7c5d48', '#735743', '#674e3d', '#5b4636', '#503d30', '#45352a', '#3a2d24', '#30261e', '#261e19', '#1c1713', '#120e0b', '#020101'] #orange sequential
#colours = ['#a8c2ff', '#9cb4ec', '#90a5d8', '#8497c6', '#7a8bb6', '#7080a6', '#6a789c', '#637192', '#5d6a89', '#55617d', '#4e5871', '#464f65', '#40485c', '#3a4153', '#353b4a', '#2f3441', '#2a2e39', '#242731', '#1f2129', '#1a1b21', '#141619', '#0d0e11', '#040405'] #blue sequential
# colours = ['#a8c2ff', '#6a789c', '#40485c', '#1a1b21', '#0d0e11', '#040405'] #blue sequential cut down
# colours = ['#a85f5f', '#a56f75', '#9f7e8a', '#99899b', '#9292a9', '#8a9bb7', '#7da6c8', '#62b5e1'] # black figure test
# colours = ['#a85f5f', '#9f7e8a', '#9292a9', '#7da6c8', '#62b5e1'] # black figure test cut down
# colours = ['#ffffff', '#e8eaf4', '#c4c9e2', '#95a1cb', '#7587bc', '#556fad', '#3b5fa3', '#104e98'] # sum black figure cut down
colours = ['#ffffff', '#fbf9ff', '#f7f3fe', '#f1ebfe', '#ebe2fe', '#e5dafd', '#ddd0fc', '#d5c6fc', '#cdbbfb', '#c5b1fa', '#bda7f9', '#b49df8', '#ac93f7', '#a28af6', '#9980f5', '#8f76f4', '#856df3', '#7963f2', '#6d59f0', '#6050ef', '#5146ee', '#3e3dec'] # sum black figure
colourIndex = 0

bg_color = 'black'
fg_color = 'white'
f1 = plt.figure(facecolor=bg_color, edgecolor=fg_color)
axes = plt.axes((0.1, 0.1, 0.8, 0.8), axisbg=bg_color)
axes.xaxis.set_tick_params(color=fg_color, labelcolor=fg_color)
axes.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color)
for spine in axes.spines.values():
    spine.set_color(fg_color)
# plt.title("Cib Power Spectra by shell", color=fg_color)                                
plt.title("Cib Power Spectra summed", color=fg_color)                                
plt.xlabel(r'$l$', fontsize=40, color=fg_color)
plt.ylabel(r'$C_l$', fontsize=40, color=fg_color)

residualCls = []
residualEll = np.arange(0, lmax+1)

while i<n:
    print "i:", i
    # if i%6 != 0:
    #     i += 1
    #     continue

    zMin = i*dz
    zMax = (i+1)*dz
    zMinCib = (i+1)*dz
    zMaxCib = (i+2)*dz

    labelBase = "%.2f<z<%.2f"%(zMinCib, zMaxCib)
    unlensedMapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_unlensed/cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_ns2048_tot.fits" % (zMinCib, zMaxCib)
    unlensedClFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_unlensed_cl/cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_ns2048_tot.dat" % (zMinCib, zMaxCib)
    unlensedEll = np.arange(0, lmax + 1)
    if os.path.exists(unlensedClFilename):
        unlensedCl = hp.read_cl(unlensedClFilename)
    else:
        unlensedMap = hp.read_map(unlensedMapFilename)
        unlensedMap = np.nan_to_num(unlensedMap)
        unlensedCl = hp.anafast(unlensedMap)
        unlensedCl[0] = 0
        unlensedCl[1] = 0
        hp.write_cl(unlensedClFilename, unlensedCl)
    unlensedCl = savitzky_golay(unlensedCl, 75, 3)
    summedCl += unlensedCl
    # if i%dPlot == 0:
        # plt.semilogy(unlensedEll[1:], unlensedCl[1:], colours[colourIndex], ls='--', label=labelBase + " unlensed", linewidth=2, axes=axes)
        # plt.semilogy(unlensedEll[1:], summedCl[1:], colours[colourIndex], ls='--', label="z<%.2f"%zMaxCib, linewidth=2, axes=axes)
    
    lensedMapFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_lensed/lensed_cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_ns2048_tot.fits" % (zMinCib, zMaxCib)
    lensedClFilename = "/scratch2/r/rbond/phamloui/lenspix_files/cib_v2_lensed_cl/lensed_cib_fullsky_ns2048_zmin%.2f_zmax%.2f_nu217_ns2048_tot.dat" % (zMinCib, zMaxCib)
    lensedEll = np.arange(0, lmax + 1)
    if os.path.exists(lensedClFilename):
        lensedCl = hp.read_cl(lensedClFilename)
    else:
        lensedMap = hp.read_map(lensedMapFilename)
        lensedMap = np.nan_to_num(lensedMap)
        lensedCl = hp.anafast(lensedMap)
        lensedCl[0] = 0
        lensedCl[1] = 0
        hp.write_cl(lensedClFilename, lensedCl)
    lensedCl = savitzky_golay(lensedCl, 75, 3)
    summedLensedCl += lensedCl
    if i%dPlot == 0:
        # plt.semilogy(lensedEll[1:], lensedCl[1:], colours[colourIndex], label=labelBase + " lensed", linewidth=2, axes=axes)
        plt.semilogy(lensedEll[1:], summedLensedCl[1:], colours[colourIndex], label="z<%.2f"%zMaxCib, linewidth=2, axes=axes)
        # colourIndex += 1
        residualCl = (lensedCl - unlensedCl) / unlensedCl
        residualCls.append((labelBase, residualCl))

    i += 1
    colourIndex += 1  

legend = plt.legend(numpoints=1, loc="lower right", shadow=True)
frame = legend.get_frame()
frame.set_facecolor(bg_color)
frame.set_edgecolor(fg_color)
for text in legend.get_texts():
    text.set_color(fg_color)

f2 = plt.figure(facecolor=bg_color, edgecolor=fg_color)
axes = plt.axes((0.1, 0.1, 0.8, 0.8), axisbg=bg_color)
axes.xaxis.set_tick_params(color=fg_color, labelcolor=fg_color)
axes.yaxis.set_tick_params(color=fg_color, labelcolor=fg_color)
for spine in axes.spines.values():
    spine.set_color(fg_color)
plt.title("Shell Power Spectra Residuals", color=fg_color)                              
plt.xlabel(r'$l$', fontsize=40, color=fg_color)
plt.ylabel(r'$\Delta C_l/C_l$', fontsize=40, color=fg_color)

_colIndex = 0
for _label, residualCl in residualCls:
    plt.plot(residualEll, residualCl, colours[_colIndex], label=_label)
    _colIndex += dPlot

legend = plt.legend(numpoints=1, loc="lower right", shadow=True)
frame = legend.get_frame()
frame.set_facecolor(bg_color)
frame.set_edgecolor(fg_color)
for text in legend.get_texts():
    text.set_color(fg_color)

plt.show()

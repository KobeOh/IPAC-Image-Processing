from astropy.io import fits
from astropy.table import Table, join
import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
from PIL import Image, ImageDraw, ImageFont
from astropy.io import ascii
import tarfile
import os
import shutil
from datetime import datetime


#INPUT VARIABLES HERE
#Each input variable should be a list with as many elements as there are galaxies to be processed.
#If a specific column does not exist, use the source id index as a dummy index

#filename of object table
allTables = ['hlsp_legus_hst_acs-wfc3_ngc5474_multiband_v1_padagb-mwext-avgapcor.tab',
             'ngc4449_uvis_working_with4s.txt']
#name of galaxy
allGals = ['ngc5474',
           'ngc4449']
#Wavelengths (as strings)
allWavelengths = [['275', '336', '438', '606', '814'],
                  ['275', '336', '435', '555', '814']]
#File Names corresponding to each wavelenght's .fits file
allFitsFiles = [['hlsp_legus_hst_uvis_ngc5474-c_f275w_v1_drc.fits', 'hlsp_legus_hst_uvis_ngc5474-c_f336w_v1_drc.fits','hlsp_legus_hst_uvis_ngc5474-c_f438w_v1_drc.fits', 'hlsp_legus_hst_acs_ngc5474-c_f606w_v1_drc.fits', 'hlsp_legus_hst_acs_ngc5474-c_f814w_v1_drc.fits'],
                ['hlsp_legus_hst_uvis_ngc4449_f275w_v1_drc.fits','hlsp_legus_hst_uvis_ngc4449_f336w_v1_drc.fits','hlsp_legus_hst_acs_ngc4449_f435w_v1_drc.fits','hlsp_legus_hst_acs_ngc4449_f555w_v1_drc.fits','hlsp_legus_hst_acs_ngc4449_f814w_v1_drc.fits']]
#Column numbers of the source id, x coordinate, y coordinate, RA, and DEC
allPosCols = [['1','2','3','4','5'],
              ['1','2','3','1','1']]
#Column numbers of the final total magnitudes in each wavelength (as strings)
allMagCols = [['6', '8', '10', '12', '14'],
              ['1', '1', '1', '1', '1']]
#Column numbers of the final photometric errors in each wavelength (as strings)
allErrCols = [['7', '9', '11', '13', '15'],
              ['1', '1', '1', '1', '1']]
#Colummn numbers of the Concentration Index, number of bands/filters, visual classification, and quality (as Strings)
allMiscCols = [['16', '33', '34', '35'],
               ['1', '1', '4', '5']]
#Wavelength of Concentration Index (as a String)
allCI = ['606',
         '555']
#Whether this galaxy has the number of bands/filters to sort by
allHasNumBands = [True,
                False]
#CONSTANT VARIABLES
#number of bands detected (This should always be 4)
minNumBands = 4
#Radius of image cutouts in pixels
stamp_radius = 149
#Whether errors are printed to a separate file
printErrors = False



def createMEFs(tableFile, gal, cat, qual, wavelengths, fitsFileNames, posCols, magCols, errCols, miscCols, CI, hasNumBands):

    t = Table.read(tableFile, format='ascii')

#Edge Case Handling for ngc5457c, which resets object ID numbers after 9999
    if gal == 'ngc5457c':
        for i in range(len(t[9999::])):
            t[9999+i]['col' + posCols[0]] = 10000+i

#sort by cluster class (cat), number of bands detected (num-band), and quality (qual)
    cut = t[(t['col' + miscCols[2]]==cat)]
    cut = cut[(cut['col' + miscCols[3]] == qual)]
#TODO: determine what a minNumBands of 200 means
    if hasNumBands:
            cut = cut[(cut['col' + miscCols[1]]>=minNumBands)]
    print(gal + ' contains ' + str(len(cut)) + ' ' + ' class ' + str(cat) + ' quality ' + str(qual) + " objects.")

#Load full images for each band

    im = []
    for i in range(len(wavelengths)):
        im.append([fitsFileNames[i],'f' + wavelengths[i] + 'w'])

#Create .fits cutouts 299x299pix centered around each object
#NOTE: .fits files store data as an array[y][]
    for i in range(len(cut)):
        for j in range(len(im)):
            obj, x, y = cut['col' + posCols[0]][i], cut['col' + posCols[1]][i], cut['col' + posCols[2]][i]
            obj, x, y = int(round(obj)), int(round(x)), int(round(y))
            image = fits.open(str(im[j][0]))
            image_data = image[0].data
            new_array = image_data[(y-stamp_radius-1):(y+stamp_radius),(x-stamp_radius-1):(x+stamp_radius)]
            new_image = fits.PrimaryHDU(new_array)
            hdulist = fits.HDUList([new_image])
            name = str(gal) + '_' + im[j][1] + '_obj_' + str(obj) + '_class' + str(cat) + '_quality' + str(qual)
            hdulist.writeto(name + '.fits', overwrite = True)

#set up MEF-format HDUs for all bands
    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU())
    new_hdul[0].name = gal

    for i in range(len(wavelengths)):
        new_hdul.append(fits.ImageHDU())
        new_hdul[1 + i].name = 'f' + wavelengths[i] + 'w'

#copy .fits data from each band into an MEF for all objects
    objects = cut['col' + posCols[0]]
    names = []
    errors = []

    for i in range(len(objects)):
        hdr = new_hdul[0].header
        hdr['Galaxy']=gal
        hdr['ObjectID']=str(int(round(objects[i])))
        hdr['Class']=str(cat)
        hdr['Quality']=str(qual)
        hdr['x_coord']=cut['col' + posCols[1]][i]
        hdr['y_coord']=cut['col' + posCols[2]][i]
        hdr['RA']=cut['col' + posCols[3]][i]
        hdr['Dec']=cut['col' + posCols[4]][i]
        for w in range(len(wavelengths)):
            hdr['m_f' + wavelengths[w] + 'w']=cut['col' + magCols[w]][i]
            hdr['e_f' + wavelengths[w] + 'w']=cut['col' + errCols[w]][i]
        hdr['CI_'+ CI]=cut['col' + miscCols[0]][i]
        hdr['N_filt']=cut['col' + miscCols[1]][i]

        for w in range(len(wavelengths)):
            wImage = fits.open(gal + '_f' + wavelengths[w] + 'w_obj_' + str(int(round(objects[i]))) + '_class' + str(cat) + '_quality' + str(qual) + '.fits')
            new_hdul[w + 1].data = wImage[0].data
            if np.min(wImage[0].data)==np.max(wImage[0].data):
                 errors.append((str(int(round(objects[i]))),'f' + wavelengths[w] +   'w'))
                 #print('error on ' + str(int(round(objects[i]))) + ', f' + wavelengths[w] + 'w = %.2e' % np.max(wImage[0].data))
        #if np.max(wImage[0].data)>=500:
        #    highpix.append((str(int(round(objects[i]))),'f' + wavelengths[w] + 'w',np.max(wImage[0].data)))
            if (cut['col' + magCols[w]][i] == 66.666):
                print(cut['col' + posCols[0]][i], ' has f' + wavelengths[w] + ' zeroed')
                for j in range(len(new_hdul[1 + w].data)):
                    for k in range(len(new_hdul[1 + w].data)):
                        new_hdul[w + 1].data[j][k] = 0.
            wImage.close()



        #Save MEF & copy names for .tar use

        name = 'MEF_' + gal + '_obj_' + str(int(round(objects[i]))) + '_class' + str(cat) + '_quality' + str(qual) +'.fits'
        names.append(name)
        new_hdul.writeto(name, overwrite = True)
        new_hdul.close()

    #remove individual .fits files after each MEF is saved
        for w in range(len(wavelengths)):
            os.remove(gal + '_f' + wavelengths[w] + 'w_obj_' + str(int(round(objects[i]))) + '_class' + str(cat) + '_quality' + str(qual) + '.fits')

    #Put MEFs into a new directory
    #TODO: Add Compression

    if errors != [] and printErrors:
        errors_table = Table(rows=errors, names=('ID', 'Filter'))
        ascii.write([errors_table['ID'], errors_table['Filter']],gal + '_class' + str(cat) + '_quality' + str(qual) + '_errors.tab', names =['ID','Filter'], overwrite='True')

    cwd = os.getcwd()
    galDir = '\\' + gal
    newDir = '\\MEF ' + gal + ' Class ' + str(cat) + ' Quality ' + str(qual)

    if cat == 1 and qual == 0:
        shutil.rmtree(cwd + galDir, True)
        os.mkdir(cwd + galDir)

    os.mkdir(cwd + galDir + newDir)


    for name in names:
        os.rename(cwd + '\\' + name, cwd + galDir + newDir + '\\' + name)

def main():
    for i in range(0,len(allTables)):
        for c in range(1,5):
            for q in range(0,3):
                createMEFs(allTables[i], allGals[i], c, q, allWavelengths[i], allFitsFiles[i], allPosCols[i], allMagCols[i], allErrCols[i], allMiscCols[i], allCI[i], allHasNumBands[i])

main()

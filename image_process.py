''' ############################ NIRC2 Image Processing ##############################
                             written by Logan Pearce, 2019
    ##################################################################################
    Prepare raw NIRC2 images for science by linearizing, dark subtracting, flat dividing,
    and interpolating over bad pixels.  Python 2.7 based. This script must be located in
    the same directory as the files nirc2.512.512.badpix and nirc2.1024.1024.badpix

# Requires:
#   python packages astropy
#   files nirc2.512.512.badpix and nirc2.1024.1024.badpix
#
# Input:
#   Directory containing NIRC2 science images
#   Directory containing associated NIRC2 dark and flat calibration images
#
# Output:
#      science_image_original_filename_LDFC: linearized, dark subtracted, flat corrected .fits file
#      science_image_original_filename_LDFBC: same as above with bad pixels interpolated over .fits file
#
# Useage:
# image_process.py [-h] [path_to_science_image_directory] [-d ANG_DIAM] [-c CONDITIONS]
#                      [-t TABLE]
#                      filename
#
# positional arguments:
# images                the path to the directory containing science images,
#                       including final '/'
# cals                  the path to directory containing calibration images,
                        including final '/'

# optional arguments:
# -h, --help            show this help message and exit
# -f MASTER_FLAT, --master_flat MASTER_FLAT
#                       input a previously made master flat image.
# -d MASTER_DARK, --master_dark MASTER_DARK
#                       input a previously made master dark image.
#
# examples:
#    python image_process.py /DSTuc/science/2017/ /DSTuc/cal/2017   <- Calibrate all science images of DS Tuc in the 2017 observation
#                            by looking through all of the calibration files from that observation and finding the right ones to apply
#    python image_process.py /DSTuc/science/2017/ /DSTuc/cal/2017 -f flat.20170702.4000.3.5.2.fits   <- Calibrate the
#                            science images by finding the appropriate dark frames, but using this previously made master flat frame.

'''

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.io.fits import getheader
import os
import time
import argparse
import warnings
warnings.filterwarnings("ignore")


def mode(a):
    '''Counts the number of unique occurances of values and returns which value occurs the most.  
    Made this because the scipy mode function was so very slow on these files.'''
    nums,counts = np.unique(a, return_counts=True)
    return nums[np.argmax(counts)]

def make_master_dark(list):
    ''' Makes a median-combined and linearized master dark file from a stack of calibration images
        Inputs:
            list (array): array of calibration filenames to put into stack
        Returns:
            masterdark (2d array): 2d image array of master dark
    '''
    # Linearizing coefficients:
    coeff = [1.001,-6.9e-6,-0.70e-10]
    # 10 darks is more than enough to make a master:
    if len(list) > 10:
        length = 10
    else:
        length = len(list)
    # Initialize master dark image:
    darkstack = np.zeros((length,imsize1,imsize2))
    # For each dark image;
    for line,i in zip(list,range(length)):
        # Open the header and pull out coadds:
        drkhdr = getheader(line, ignore_missing_end=True)
        coadds = float(drkhdr['coadds'])
        # Open the data file:
        dark1i = fits.open(line, ignore_missing_end=True)
        dark1s = dark1i[0].data
        # Linearize the dark frame:
        norm1 = coeff[0]+(coeff[1]*(dark1s/coadds))+(coeff[2]*((dark1s/coadds)**2))
        dark1 = dark1s/norm1
        # Place into the stack:
        darkstack[i]= dark1

    # Median combine the corrected dark frames into master dark frame:
    masterdark = np.median(darkstack,axis=0)
    fits.writeto('masterdark.fits',masterdark,drkhdr,overwrite=True)
    return masterdark
    
def make_master_flat(list,masterdark):
    ''' Makes a median-combined, linearized, and dark subtracted master flat from a stack of calibration images
        Inputs:
            list (array): array of calibration filenames to put into stack
        Returns:
            masterflat (2d array): 2d image array of master flat
    '''
    # Linearizing coefficients:
    coeff = [1.001,-6.9e-6,-0.70e-10]
    # 10 flats is more than enough to make a master:
    if len(list) > 10:
        length = 10
    else:
        length = len(list)
    # Initialize flat stack:
    flatstack = np.zeros((length,imsize1,imsize2))
    # For each flat frame:
    for line,i in zip(list,range(length)):
        # Open the header and get the coadds:
        flthdr = getheader(line, ignore_missing_end=True)
        coadds = float(flthdr['coadds'])
        # Open the image:
        flat1i = fits.open(line, ignore_missing_end=True)
        flat1s = flat1i[0].data
        # Linearize:
        norm1 = coeff[0]+(coeff[1]*(flat1s/coadds))+(coeff[2]*((flat1s/coadds)**2))
        flat1l = flat1s/norm1
        # Dark subtract the flat frame:
        flat1=flat1l-masterdark
        # Normalize the flat:
        flat1 = flat1/mode(flat1)
        # Add it to the stack:
        flatstack[i]= flat1
    # Median combine the stack:
    masterflat = np.median(flatstack,axis=0)
    # Normalize the master:
    masterflat = masterflat / np.mean(masterflat)
    # Write out calibrated science frame:
    fits.writeto('masterflat.fits',masterflat,flthdr,overwrite=True)
    return masterflat

def calibrate_science_image(image,master_dark,master_flat,newfilename):
    ''' Makes linearizes, dark-subtracts, and flat-divides science frames
        Inputs:
            image (2d array): image to be calibrated
            masterdark,masterflat (2d array): master dark and flat arrays
            newfilename (str): name for calibrated image file
        Returns:
            flat_div (2d array): calibrated image
            imhdr (fits header): annotated header
    '''
    # Linearizing coefficients:
    coeff = [1.001,-6.9e-6,-0.70e-10]
    # Open science image:
    im1 = fits.open(image, ignore_missing_end=True)
    im = im1[0].data
    # Get coadds:
    imhdr = getheader(image, ignore_missing_end=True)
    coadds = float(imhdr['coadds'])
    # Linearize science frame:
    norm1 = coeff[0]+(coeff[1]*(im/coadds))+(coeff[2]*((im/coadds)**2))
    im = im/norm1
    # Dark subtract and flat divide image:
    flat_div = (im - master_dark) / master_flat
    # Add a comment to the header:
    imhdr['COMMENT'] = '         Linearized, dark subtracted and flat corrected on '+time.strftime("%m/%d/%Y")+ ' By Logan A. Pearce'
    # Write out calibrated science frame:
    fits.writeto(newfilename,flat_div,imhdr,overwrite=True)
    return flat_div, imhdr

def bad_pixel_fix(image,hdr,sizex,sizey,newfilename,s=3):
    ''' Fixes bad pixels by interpolating from surrounding pixels.
        Copied from:
        http://mtham.ucolick.org/egates/2016GradWorkshop/PDFs/DataReduction/DataReductionProcedures-Python-2016.pdf
        Inputs:
            image (2d array): image to be fixed
            hdr (fits header): header of image for getting
            sizex, sizey (float): dimensions of image to be fixed
            newfilename (str): new file name
            s (int): radius of pixels to use in interpolation (default = 3)
        Returns:
            flat_div (2d array): pixel corrected image
    '''

    # Depending on the size of the image, load the correct bad pixel list:
    #     (Python indicies are (row,column) so the image indicies are (y,x) in python; additionally python begins with
    #      initial index of 0, while image data begins at 1)
    if sizex == 1024 and sizey == 1024:
        badpix = np.loadtxt(open("nirc2.1024.1024.badpix","rb"))
        badpix = badpix.astype(int)
    elif sizex == 512 and sizey == 512:
        badpix = np.loadtxt(open("nirc2.512.512.badpix","rb"))
        badpix = badpix.astype(int)
    else:
        print "I don't have a bad pixel list to match this image size.  I only have 1024x1024 or 512x512.  \
        Check your image sizes please."
        quit()

    # Make a boolean mask the same size as the image and set all initial values to False:
    mask = np.ma.make_mask(image,copy=True,shrink=True, dtype=np.bool)
    mask[:,:] = False
    # For each bad pixel in the list, set the value of the mask to true:
    for i in badpix:
        mask[i[1]-1][i[0]-1] = True
    # Set the value of all bad pixels in the image to nan:
    mdata = np.ma.masked_array(image,mask=mask,fill_value=np.nan)
    # Make a new array as a copy of original image:
    badpixelfixed = image.copy()
    # For each x and y value, loop through image and replace the value of all "nan"
    # pixels with the mean of the four pixels on either side of the bad one:
    for i in range(0,mdata.shape[0]):
        for j in range(0,mdata.shape[1]):
            if np.math.isnan(mdata[i,j]):
                x1 = i-s
                x2 = i+s+1
                y1 = j-s
                y2 = j+s+1
                if x1<0:
                    x1 = 0
                if x2>mdata.shape[0]:
                    x2=mdata.shape[0]
                if y1<0:
                    y1 = 0
                if y2>mdata.shape[1]:
                    y2 = mdata.shape[1]
                badpixelfixed[i,j] = np.mean(mdata[x1:x2,y1:y2])
    # Add a comment to the header:
    hdr['COMMENT'] = '         Bad pixels fixed on '+time.strftime("%m/%d/%Y")+ ' By Logan A. Pearce'
    # Write out the bad pixel corrected image to a new fits file:
    fits.writeto(newfilename,badpixelfixed,hdr,overwrite=True)
    return badpixelfixed
    
    


# Pull out arguments:
parser = argparse.ArgumentParser()
# Required positional arguments:
parser.add_argument("images", help="the path to the directory containing science images, including final '/' ", type=str)
parser.add_argument("cals", help="the path to directory containing calibration images, including final '/' ", type=str)
# Optional positional arguments"
parser.add_argument("-f","--master_flat", help="input a previously made master flat image.",type=str)
parser.add_argument("-d","--master_dark", help="input a previously made master dark image.",type=str)

args = parser.parse_args()
# Point at the directory containing the images to processed:
images_directory = args.images
# Point at the directory containing all calibration files:
cals_directory = args.cals

os.system('ls '+images_directory+'*.fits > images_list')
os.system('ls '+cals_directory+'*.fits > cals_list')

with open('images_list') as f:
    images = f.read().splitlines()
with open('cals_list') as f:
    cals = f.read().splitlines()


# For each science image:
for image in images:
    print 'Calibrating ',image
    # Open the header:
    imhdr = getheader(image, ignore_missing_end=True)
    # Pull out the relevant header parameters:
    mm,ss,ii,cc= imhdr['MULTISAM'],imhdr['SAMPMODE'],imhdr['ITIME'],imhdr['COADDS']
    imfilt=imhdr['FILTER']
    #print imfilt
    imsize1,imsize2 = imhdr['NAXIS1'],imhdr['NAXIS2']
    #print imsize1,imsize2
    try:
        imkoaid = imhdr['KOAID']
        obsdate = imkoaid.split('.')[1]
    except:
        imkoaid = imhdr['OBJECT'].replace(" ","_")
        obsdate = imhdr['DATE-OBS'].replace("-", "")
    print imkoaid,obsdate
    

    if args.master_dark:
        pass
    else:
        # Go through each calibration file and find the matching cal files:
        darklist = []
        flatlist = []
    
        for cal in cals:
            hdr = getheader(cal, ignore_missing_end=True)
            m,s,itime,c = hdr['MULTISAM'],hdr['SAMPMODE'],hdr['ITIME'],hdr['COADDS']
            filt = hdr['FILTER']
            print filt
            try:
                caltype = hdr['KOAIMTYP']
                koaid = hdr['KOAID']
                caldate = koaid.split('.')[1]
            except:
                caltype = hdr['OBJECT']
                caldate = hdr['DATE-OBS'].replace("-", "")
            size1,size2 = hdr['NAXIS1'], hdr['NAXIS2']
            #print size1,size2
            if 'dark' in caltype:
                # If dark, match multisam, sampmode, itime, coadds, and the date
                if m==mm and s==ss and itime==ii and c==cc and size1==imsize1 and size2==imsize2 and caldate==obsdate:
                    darklist.append(cal)
            if 'flat' in caltype:
                # If flat, match filter and date:
                if filt==imfilt and size1==imsize1 and size2==imsize2: #and caldate==obsdate:
                    flatlist.append(cal)
    print darklist
    print flatlist
    # Make master dark:
    if args.master_dark:
        masterdark = fits.open(args.master_dark)[0].data
    else:    
        masterdark = make_master_dark(darklist)

    # Make master flat:
    if args.master_flat:
        masterflat = fits.open(args.master_flat)[0].data
    else:
        masterflat = make_master_flat(flatlist,masterdark)
    
    # Calibrate science frame:
    newfilename = '..'+image.split('.')[2]+'.'+image.split('.')[3]+'.'+image.split('.')[4]+'.LDFC.'+image.split('.')[-1]
    image_ldfc, image_ldfc_hdr = calibrate_science_image(image,masterdark,masterflat,newfilename)

    # Fix bad pixels:
    newfilename = '..'+image.split('.')[2]+'.'+image.split('.')[3]+'.'+image.split('.')[4]+'.LDFBC.'+image.split('.')[-1]
    image_ldfbc = bad_pixel_fix(image_ldfc,image_ldfc_hdr,size1,size2,newfilename,s=3)
    print 'done... next'

os.system('rm images_list')
if args.master_dark:
    pass
else:
    os.system('rm cals_list')

    
    


            





###############################################################################
#                                                                             #
#                        NIRC2 Image Processing pipeline                      #
#                                                                             #
#                        Written by Logan A. Pearce (2017)                    #
###############################################################################

# Search through a folder of calibration images and find the ones that match your science images.
# Writes out a list of the file names for the dark frames and flat fields frames that match
# the science image.
# Input:
#   - Path + filename of science image to match the dark and flat frames to
#   - Path to folder containing calibration images to search for matches
# Output:
#   - List of dark images that match the science image's Multisam, Sampmode, ITime, and Coadds called "darklist"
#   - List of flat images that match the science image's filter called "flatlist"

# execution syntax:
# python finding_darkflat.py path_to_science_image path_to_calibration_images (do not include final "/" in path)
# example:
# python finding_darkflat.py KOA_11708_sci/NIRC2/raw/sci/N2.20170628.26610.fits KOA_11708_cal/NIRC2/raw/cal
# 
# First: make a list of all .fits files in the folder containing calibration images:
#     ls ~/Desktop/UTexas/Astro_research_data/data/GSC6214/2017_06_27/KOA_23637_cal/NIRC2/raw/cal/*.fits>list
# Next: execute "finding_darkflat.py" to search through those files for the ones that match the science frame:
#     python finding_darkflat.py ~/Desktop/UTexas/Astro_research_data/data/GSC6214/2017_06_27/N2.20170628.26771.fits
# Repeat for all calibration folders.  This script outputs two files in the science fram folder titled "darklist" and "flatlist"
# Next: Make a list of all the science frame images titled "objectlist" and put into the folder containing science frames:
#     ls ~/Desktop/UTexas/Astro_research_data/data/GSC6214/2017_06_27/*.fits>~/Desktop/UTexas/Astro_research_data/data/GSC6214/2017_06_27/objectlist
# Last: run the image process script on all science frames:
#     python image_process.py ~/Desktop/UTexas/Astro_research_data/data/GSC6214/2017_06_27/


from astropy.io import fits
from astropy.io.fits import getheader
import os

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("image_filename",type=str)
#parser.add_argument("path",type=str)
args = parser.parse_args()
#path=args.path
image=args.image_filename

#folder='~/'+path.split('/')[2]+'/'+path.split('/')[3]+'/'+path.split('/')[4]+'/'+path.split('/')[5]+'/'+path.split('/')[6]+'/'+path.split('/')[7]+'/' \
#  +path.split('/')[8]+'/' +path.split('/')[9]+'/' +path.split('/')[10]+'/' +path.split('/')[11]+'/'
directory = image.split('/')[0]+'/'+image.split('/')[1]+'/'+image.split('/')[2]+'/'+image.split('/')[3]+'/'+image.split('/')[4]+'/'+image.split('/')[5]+ \
  '/'+image.split('/')[6]+'/'+image.split('/')[7]+'/'+image.split('/')[8]+'/'

#os.system('ls '+folder+'*.fits > list')

# Make a temporary list of all calibration images in the specified folder:
#os.system('ls '+path+'*.fits > list')

# Open the science image to find the parameters needed to match in the darks and flats:
print "Opening image ", args.image_filename
imhdr = getheader(args.image_filename)
mm,ss,ii,cc= imhdr['MULTISAM'],imhdr['SAMPMODE'],imhdr['ITIME'],imhdr['COADDS']
imfilt=imhdr['FILTER']
imsize = imhdr['NAXIS1']
print 'Searching for:'
print 'Multisam: ',mm,'Sampmode: ',ss,'ITime: ',ii,'Coadds: ',cc
print 'Filter: ',imfilt
print 'Image size: ',imsize

with open('list') as f:
    z = f.read().splitlines()
a = open(directory+'darklist', 'a')
b = open(directory+'flatlist', 'a')

for line in z:
    hdr = getheader(line)
    m,s,i,c = hdr['MULTISAM'],hdr['SAMPMODE'],hdr['ITIME'],hdr['COADDS']
    obj,filt=hdr['object'],hdr['FILTER']
    imtype = hdr['KOAIMTYP']
    size = hdr['NAXIS1']
    #print 'obj',obj
    # In the calibration image is a dark frame, match the multisam, samp mode, integration time, and number
    # of coadds:
    if obj=='darks' or obj=='really_dark' or obj=='dark' or obj=='Dark':
        print m,s,i,c
        print hdr['NAXIS1']
        if m==mm and s==ss and i==ii and c==cc and size==imsize:
             a.write(str(line) + "\n")
        else:
            pass
    if imtype=='darks' or obj=='really_dark' or obj=='dark' or obj=='Dark':
        print m,s,i,c
        print hdr['NAXIS1']
        if m==mm and s==ss and i==ii and c==cc and size==imsize:
             a.write(str(line) + "\n")
        else:
            pass
    # if the calibration image is a flat field, match the filter type:
    elif obj=='LampOn' or obj=='LampOff' or obj=='really_dome' or obj=='dome' or obj=='dome flat' or imtype=='flatlamp':
        print filt
        print hdr['NAXIS1']
        if filt==imfilt and size==imsize:
            b.write(str(line) + "\n")
        else:
            pass
    else:
        pass

a.close()
b.close()
print 'writing dark and flat list... done'

#os.system('rm list')


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
# Written by Logan A. Pearce

from astropy.io import fits
from astropy.io.fits import getheader
import os

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("image_filename",type=str)
parser.add_argument("path",type=str)
args = parser.parse_args()
path = args.path

# Make a temporary list of all calibration images in the specified folder:
os.system('ls '+path+'/*.fits > list')

# Open the science image to find the parameters needed to match in the darks and flats:
print "Opening image ", args.image_filename
imhdr = getheader(args.image_filename)
mm,ss,ii,cc= imhdr['MULTISAM'],imhdr['SAMPMODE'],imhdr['ITIME'],imhdr['COADDS']
imfilt=imhdr['FILTER']
print 'Searching for:'
print 'Multisam: ',mm,'Sampmode: ',ss,'ITime: ',ii,'Coadds: ',cc
print 'Filter: ',imfilt

with open('list') as f:
    z = f.read().splitlines()
a = open('darklist', 'w')
b = open('flatlist', 'w')

for line in z:
    hdr = getheader(line)
    m,s,i,c = hdr['MULTISAM'],hdr['SAMPMODE'],hdr['ITIME'],hdr['COADDS']
    obj,filt=hdr['object'],hdr['FILTER']
    #print 'obj',obj
    # In the calibration image is a dark frame, match the multisam, samp mode, integration time, and number
    # of coadds:
    if obj=='darks' or obj=='really_dark' or obj=='dark' or obj=='Dark':
        if m==mm and s==ss and i==ii and c==cc:
             a.write(str(line) + "\n")
        else:
            pass
    # if the calibration image is a flat field, match the filter type:
    elif obj=='LampOn' or obj=='LampOff' or obj=='really_dome' or obj=='dome' or obj=='dome flat':
        #print 'filt',filt
        if filt==imfilt:
            b.write(str(line) + "\n")
        else:
            pass
    else:
        pass

a.close()
print 'writing dark and flat list... done'

os.system('rm list')


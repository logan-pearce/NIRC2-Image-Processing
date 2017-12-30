###############################################################################
#                                                                             #
#                        NIRC2 Image Processing pipeline                      #
#                                                                             #
#                        Written by Logan A. Pearce (2017)                    #
###############################################################################

# Search through a folder of science images for the images of the specific object you are looking for.
# Writes out a list of the filenames for the images that list that object in the header.
# Input: object name, folder containing science images
# Output:
#    - list of science images of the desired object called "objectlist"
#    - list of each image's parameters for finding the corresponding dark and flat frames for processing
#        called "masterlist"

# execution syntax:
# python finding_object.py object path_to_science_images (do not include final "/" in path)
# example:
# python finding_object.py GSC6214-210 KOA_11708_sci/NIRC2/raw/sci
#
# Written by Logan A. Pearce

from astropy.io import fits
from astropy.io.fits import getheader
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("object",type=str)
parser.add_argument("path",type=str)
args = parser.parse_args()
path = args.path

# Create a new file to output filelist to:
a = open('objectlist', 'w')

# Create a temporary list of all science fits images in the folder:
os.system('ls '+path+'/*.fits > scilist')

#Open the list and read each line:
print 'Searching science images for desired object...'
with open('scilist') as f:
    z = f.read().splitlines()

# For each image in the list, look for ones with the specified object in the header:
for line in z:
    hdr = getheader(line)
    obj=hdr['object']
    if obj==args.object:
        a.write(str(line) + "\n")

a.close()
print 'Found desired image files... done'

# Delete temporary file:
os.system('rm scilist')

# Generate of all found science image files and their relevant parameters for matching dark frames and flat fields:
print 'Opening science images and creating list of dark and flat parameters...'

with open('objectlist') as f:
    z = f.read().splitlines()

for line in z:
    hdr = getheader(line)
    imagename = line.split('.')[2]
    m,s,i,c = hdr['MULTISAM'],hdr['SAMPMODE'],hdr['ITIME'],hdr['COADDS']
    obj,filt=hdr['object'],hdr['FILTER']
    string = imagename + ' Multisam: '+str(m)+' Sampmode: '+str(s)+' ITime: '+str(i)+' Coadds: '+str(c)+' Filter: '+str(filt)+' '
    k = open('masterlist', 'a')
    k.write(string + "\n")
    k.close()

print 'All done.'

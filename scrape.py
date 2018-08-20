###############################################################################
#                                                                             #
#                        NIRC2 Image Processing pipeline                      #
#                                  scrape.py                                  #
#                                                                             #
#                        Written by Logan A. Pearce (2017)                    #
###############################################################################
#
# Retrieve NIRC2 image files from KOA archive
# Adapted from scrape.pro by Adam Kraus.
#
# Requires: astropy, wget
#
# Input:ascii table download from KOA archive search results
#
# Output: downloaded NIRC2 fits files from KOA archive
#
# usage: scrape.py [-h] [-o OUTPUT_PATH] fitslist

#optional arguments:
#  -h, --help            show this help message and exit
#  -o OUTPUT_PATH, --output_path OUTPUT_PATH
#                        optional specify output path. Otherwise files will be
#                        downloaded to the same directory as this script.

from astropy.io import ascii
import sys
import wget
import os

############### Define arguments:
import argparse
parser = argparse.ArgumentParser()
# Required positional arguments:
parser.add_argument("fitslist", help="ascii table downloaded from KOA archive search results")
# Optional positional arguments"
parser.add_argument("-o","--output_path", help="optional specify output path.  Otherwise files will be downloaded to the same directory \
    as this script.")
args = parser.parse_args()

if args.output_path:
    out = output_path
    if os.path.exists(out):
        pass
    else:
        os.sys('mkdir '+str(out))
else:
    out=''

filein=args.fitslist
data = ascii.read(filein)
for i in range(len(data)):
    print('')
    print('Getting file ',i+1,' of ',len(data))
    query='https://koa.ipac.caltech.edu/cgi-bin/getKOA/nph-getKOA?filehand='+data['filehand'][i]+'&instrument=N2'
    if os.path.exists(str(out)+data['koaid'][i]):
        print("file exists")
    else:
        wget.download(query, str(out)+data['koaid'][i])
print('')
print('done.')

    


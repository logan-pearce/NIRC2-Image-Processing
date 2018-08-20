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
# Requires: astropy, pandas, wget
#
# Input: either ascii table download from KOA archive search results or text file
#      list of images to download
#
# Output: downloaded NIRC2 fits files from KOA archive

#
# usage: scrape.py [-h] [-a ASCII_FITSLIST] [-t TEXT_FITSLIST] [-o OUTPUT_PATH]

#optional arguments:
#  -h, --help            show this help message and exit
#  -a ASCII_FITSLIST, --ascii_fitslist ASCII_FITSLIST
#                        optional input an ascii table downloaded from KOA
#                        search results.
#  -t TEXT_FITSLIST, --text_fitslist TEXT_FITSLIST
#                        optional input a plain text file of just a list of
#                        NIRC2 filenames to be downloaded.
#  -o OUTPUT_PATH, --output_path OUTPUT_PATH
#                        optional specify output path. Otherwise files will be
#                        downloaded to the same directory as this script.

from astropy.io import ascii
import pandas as pd
import sys
import wget
import os

############### Define arguments:
import argparse
parser = argparse.ArgumentParser()
# Required positional arguments:
#parser.add_argument("fitslist", help="ascii table downloaded from KOA archive search results")
# Optional positional arguments"
parser.add_argument("-a","--ascii_fitslist", help="optional input an ascii table downloaded from KOA search results.")
parser.add_argument("-t","--text_fitslist", help="optional input a plain text file of just a list of NIRC2 filenames to be downloaded.")
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


if args.ascii_fitslist:
    filein=args.ascii_fitslist
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

elif args.text_fitslist:
    filein=args.text_fitslist
    data = pd.read_table(filein,header=None,usecols=[0],names=['koaid'])
    for i in range(len(data)):
        print('')
        print('Getting file ',i+1,' of ',len(data))
        query='https://koa.ipac.caltech.edu/cgi-bin/getKOA/nph-getKOA?filehand=/koadata9/NIRC2/'+data['koaid'][i].split('.')[1]+'/lev0/'+data['koaid'][i]+'&instrument=N2'
        if os.path.exists(str(out)+data['koaid'][i]):
            print("file exists")
        else:
            wget.download(query, str(out)+data['koaid'][i])
    print('')
    print('done.')
else:
    print('Error: Need either ascii or text fitslit input')
    parser.print_help()
    sys.exit(0)
    


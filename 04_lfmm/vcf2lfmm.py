#-*- coding: utf-8 -*-

# Input: vcf file from freebayes

# Output: matrix compatible with the lfmm format

# Help: python vcf2lfmm.py --help

#AUTHOR : LEO MILHADE (milhadel)

###################################################################
#Import libraries
import sys
import os
import numpy as np

import argparse

from scipy.stats import chisquare

import multiprocessing as mp

import warnings
warnings.filterwarnings("error")

###################################################################
# Define parser
parser = argparse.ArgumentParser(description='Converts a vcf file into a matrix compatible with the lfmm format')

parser.add_argument('input_file', type=str,
                   help='vcf file, freebayes format')

parser.add_argument('-o', '--output', dest='output_file', default="parsed_vcf.lfmm",
                    help='output file name')

args = parser.parse_args()

###################################################################
#Set parameters
###################################################################
#Set parameters
comment_char       = "#"        #comment character for header in vcf file
split_char         = "\t"       #sample field separator in vcf file
field_sep          = ":"        #for exemple in GT:DP:AD:RO:QR:AO:QA:GL field separator = ":"
DP_field           = 1          #for exemple in GT:DP:AD:RO:QR:AO:QA:GL DP field = 2
AD_field           = 2          #for exemple in GT:DP:AD:RO:QR:AO:QA:GL AD field = 3
GT_field           = 0          #for exemple in GT:DP:AD:RO:QR:AO:QA:GL AD field = 3
AD_field_sep       = ","        #separator in AD field
 

end_signal = "PROCESSING DONE"


###################################################################

def process(line):

    line = line.split(split_char)

    #Get chrom, pos, ref and alt fields

    newline = split_char.join(line[i] for i in range(0,9))

    count = []
    col_cnt = 0

    for index in range(int(indx_columns_start),nbcol) :
        col_cnt += 1
        section = line[index].split(field_sep)
        if section[GT_field] == ".":
            count += [9]
        elif section[GT_field] == "0/0":
            count += [0]
        elif section[GT_field] == "0/1":
            count += [1]
        elif section[GT_field] == "1/1":
            count += [2]
        else:
            count += [9]

    newline = split_char.join((newline, split_char.join([str(i) for i in count])))

    return newline

###################################################################


try:
    with open(args.input_file, "rU") as in_file, open(args.output_file, "w") as out_file:

        line = next(in_file)
        while line.startswith(comment_char) :
            prev_line =  line
            line = next(in_file)
        
        header = prev_line[1:]
        out_file.write(header + "\n")

        indx_columns_start = 9         #index of colum for the first sample if colum number 1 = CHROM
        nb_samples         = len(header.split(split_char)) - indx_columns_start
        nbcol              = int(indx_columns_start)+int(nb_samples)

        out_file.write(process(line) + "\n")
        for line in in_file:
            out_file.write(process(line) + "\n")

except IOError, e:
    print "File opening failed: %s" % e.strerror



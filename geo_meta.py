#!/usr/bin/env python

'''
    Generate fields needed for meta file submitted to CEO database
'''

import os
import argparse
import hashlib

parser = argparse.ArgumentParser(description='Process counts and raw data files to CEO meta table')
parser.add_argument('--instrument_model', dest = 'instrument', default = 'Illumina HiSeq 2500')
parser.add_argument('--read_length', dest = 'read_length', default = '126')
parser.add_argument('--read_layout', dest = 'read_layout', default = "paired-end")
parser.add_argument('--counts_dir', dest='counts_dir', default = '/Users/feifei/Projects/Others/Jingyi/Rcodes/counts',
                    help='specify dir where counts files are')
parser.add_argument('--rawdata_dir', dest='rawdata_dir', default = '/Volumes/molev-32-59.icm.uu.se/data/Others/Jingyi/data',
                    help='specify dir where raw data files are .fastq.gz')

# parser.add_argument('--paired_end', dest='read_layout', action='store_true')
# parser.add_argument('--single', dest='read_layout', action='store_false')
# parser.set_defaults(read_layout=True)

args = parser.parse_args()
instrument = args.instrument
read_length = args.read_length
counts_dir = args.counts_dir
rawdata_dir = args.rawdata_dir
read_layout = args.read_layout


def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

# Raw counts info
print '############ RAW COUNTS INFO ############'
for _, _, files in os.walk(counts_dir):
    for f in sorted(files):
        if f.endswith('.counts'):
            checksum = md5("/".join([counts_dir,f]))
            file_type = 'Raw read counts per gene'
            print '\t'.join([f, file_type, checksum])


# Raw data info
print '\n############### RAW DATA INFO ###########'
for _, _, files in os.walk(rawdata_dir):
    for f in sorted(files):
        if f.endswith('.fastq.gz'):
            checksum = md5("/".join([rawdata_dir,f]))
            file_type = 'fastq'               
            print '\t'.join([f, file_type, checksum, instrument, read_length, read_layout])



# Get read length, average insert size, standard deviation for paired-end experiments

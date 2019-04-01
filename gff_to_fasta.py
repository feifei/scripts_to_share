#!/usr/bin/env python
'''
    Extract fasta file from gff and genome fasta file
    No intron or whatsoever is taken care.
'''


import argparse
import sys, os, re
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from artemis_string import *

def main():
    parser = argparse.ArgumentParser(description='Extract gene fasta file from gff and genome file')
    parser.add_argument('gff_file')
    parser.add_argument('scf_file')
    parser.add_argument('--feature', dest='feature', default = "gene", type = str)
    parser.add_argument('--code', dest='code', default = 1, type = int)
    args = parser.parse_args()

    gff_file = args.gff_file
    scf_file = args.scf_file
    feature = args.feature
    code = args.code

    basename, extension = os.path.splitext(gff_file)
    nt_outfile = "%s.%s.fa" %(basename, feature)
    aa_outfile = "%s.%s.faa" %(basename, feature)


    genomes = defaultdict(Seq)
    with open(scf_file, 'r') as gh:
        for record in SeqIO.parse(gh, 'fasta'):
            genomes[record.id] = record.seq

    def get_na_seq(genomes_d, scfid, start, end, strand):
        seq = ""
        if scfid in genomes_d:
            seq = genomes_d[scfid][start-1 : end]
            seq = seq.reverse_complement() if strand == "-" else seq
        return seq
    
    def get_aa_seq(genomes_d, scfid, start, end, strand, code = 1):
        na_seq = get_na_seq(genomes_d, scfid, start, end, strand)
        return na_seq.translate(code)

    
    with open(nt_outfile, 'w') as nt_outh, open(aa_outfile, 'w') as aa_outh, open(gff_file, 'r') as gff_h:
        for line in gff_h:
            line = line.strip()
            if len(line) == 0 or re.match("#", line):
                continue
            if re.match(">", line):
                break
            arr = line.split()
            scfid = arr[0]
            type = arr[2]
            start = int(arr[3])
            end = int(arr[4])
            strand = arr[6]
            
            if type == "gene":
                geneid = re.search("ID=(.*?);", line).group(1)
                if not re.search("description=(.*?);", line):
                    description = re.search("description=(.*?)$", line).group(1)
                else:
                    description = re.search("description=(.*?);", line).group(1)
                description = filter_string(description)
            if type == feature:            
                na_seq = get_na_seq(genomes, scfid, start, end, strand)
                aa_seq = get_aa_seq(genomes, scfid, start, end, strand)
                
                print >>nt_outh, ">%s %s" %(geneid, description)
                print >>nt_outh, na_seq
                print >>aa_outh, ">%s %s" %(geneid, description)
                print >>aa_outh, aa_seq
            
if __name__ == "__main__":
    main()


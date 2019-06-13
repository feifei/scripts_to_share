'''
    Substract a list of bases to change from 
    BWA mapped Illumina data -> pileup
    - Major SNPs > 50% 
    - Indels 
'''

import re, os
import argparse
from collections import defaultdict, Counter
from Bio.SeqUtils import GC
from Bio import SeqIO
from Bio.Seq import Seq


def indel(seq):
    if not seq:
        return None
    
    seq = seq.upper()
    m = re.findall("[+-]\d+[A-Za-z]+", seq)
    return Counter(m)



def get_indel(seq, cov):
    indel_c = indel(seq)
    indel_info = []
    if indel_c:
        if len(indel_c) >= 2:
            # Check the alternative indel. In some cases, both indel has two major candidates
            (s1, c1), (s2, c2) = indel_c.most_common(2)
            if c1/float(cov) > 0.3 and cov > 10 and c2/float(cov) > 0.3:
                print scfid, pos, cov, c1/float(cov), indel_c
            
        s, c = indel_c.most_common(1)[0]
        tag = s[0]
        size = int(s[1])
        bases = s[2:]
        if c / float(cov) > 0.42 and cov >= 10:
            # 0.43 is to make 1994126 to accept +1T
            indel_info = [tag, c, bases]
    
    return indel_info
    

def base_count(seq):
    if not seq:
        return None
    
    seq = seq.upper()
    seq = re.sub("[+-]\d+[A-Za-z]+", "", seq) # Delete indels
    d = defaultdict(int)
    d["A"] = seq.count("A")
    d["C"] = seq.count("C")
    d["G"] = seq.count("G")
    d["T"] = seq.count("T")
    return d


def get_alt_base(seq, cov):
    # If base coverage >=10 X and >50%, then update to the alt_base
    count_d = base_count(seq)
    if not count_d:
        return None
    alt_base = []
    for base,c in count_d.iteritems():
        p = c/float(cov)
        if p > 0.5 and cov >= 10:
            alt_base = [base, c, cov, p]
    return alt_base
            



parser = argparse.ArgumentParser(description='Update reference seq from pileup results')
parser.add_argument('pileup_file')
parser.add_argument('--cov', dest='cov_cutoff', default=10, type=int,
                    help='Site coverage cutoff')
parser.add_argument('--percent', dest='percent', default=0.5, type=float,
                    help='Percent of alternative base cutoff')

args = parser.parse_args()
pileup_file = args.pileup_file
basename, extension = os.path.splitext(args.pileup_file)
changes_outfile = basename + ".changes.tab"

with open(changes_outfile, 'w') as outh, open(pileup_file) as pileup_inh:
    for line in pileup_inh:
        arr = line.split("\t")
        if len(arr) != 6:
            print line
            continue
        scfid, pos, ref_base, cov, seq, qual = arr
        pos = int(pos)
        cov = int(cov)
        indel_info = get_indel(seq, cov)
        if indel_info:
            tag, c, bases = indel_info
            if tag == "+":
                print >>outh, "\t".join(map(str, [scfid, pos, ".", bases, "%s:%d" %(cov, c)]))
            elif tag == "-":
                # pos+1 being the start position of deletion
                print >>outh, "\t".join(map(str, [scfid, pos + 1, bases, ".", "%s:%d" %(cov, c)])) 
        
        alt_base = get_alt_base(seq, cov)
        if alt_base:
            base, c, cov, p = alt_base
            print >>outh, "\t".join(map(str, [scfid, pos, ref_base, base, "%s:%d" %(cov, c)])) 
        

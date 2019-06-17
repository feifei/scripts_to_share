''' Parse samtools mpileup -Bf pileup results
'''
from optparse import OptionParser
import re, os
from collections import defaultdict


def clean_seq(seq):
    s = re.sub("[+-]\d+[A-Za-z]+", "", seq)
    s = re.sub("\^.|\$|\,|\.", "", s)
    return s


def get_count(seq):
    d = defaultdict(int)
    if not seq:
        return None
    d["A"] = seq.count("A") + seq.count("a")
    d["C"] = seq.count("C") + seq.count("c")
    d["G"] = seq.count("G") + seq.count("g")
    d["T"] = seq.count("T") + seq.count("t")
    
    return d

# Make it part of input 0.1 and 20
def get_alt_bases(count_d, cov):
    s = ""
    for k, v in count_d.iteritems():
        p = v/float(cov)
        if p >= 0.1 and cov >= 20:
            s += "%s:%d:%.2f;" %(k, v, p)
    return s
    

parser = OptionParser("usage: %prog mpileup_file")
(options, args) = parser.parse_args()

pileup_file = args[0]
basename, extension = os.path.splitext(pileup_file)
outfile = basename + ".snps.tab"

with open(outfile, 'w') as outh:
    with open(pileup_file, 'r') as inh:
        for line in inh:
            arr = line.split("\t")
            if len(arr) != 6:
                print line
                continue
            scfid, pos, ref_base, cov, seq, qual = arr
            c_seq = clean_seq(seq)
            count_d = get_count(c_seq)
            if not count_d:
                continue
            alt_bases = get_alt_bases(count_d, int(cov))
            if not alt_bases:
                continue
            print >>outh, "\t".join([scfid, pos, ref_base, cov, alt_bases])
        
        
        

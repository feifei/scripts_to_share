'''
    Vcf stats
    SNPs/MNPs = ?
    Insertions = ?
    Deletions = ?
    CNVs = ?
'''
import argparse
import sys, re
from collections import defaultdict
from rpy2.robjects import *

parser = argparse.ArgumentParser()
parser.add_argument('variants_vcf', type=argparse.FileType('r'), default=sys.stdin)
args = parser.parse_args()

insertions = []
deletions = []
snps = []
mnps = []
total = 0
for line in args.variants_vcf:
    if re.match("#", line):
        continue
    arr = line.split()
    ref = arr[3]
    alt = arr[4] 
    if alt == "." :
        continue
    if re.match("D", alt) or (len(alt) < len(ref)):
        # deletion
        deletions.append(ref)
    elif re.match("I", alt) or (len(alt) > len(ref)):
        # insertion
        insertions.append(alt[1:])
    else:
        # substitutions
        if len(alt) == 1:
            snps.append([alt, ref])
        else:
            mnps.append([alt, ref])
    total += 1

print "Total = %d" %total
print "Insertions = % d" %(len(insertions))
print "Deletions = % d" %(len(deletions))
print "SNPs = % d" %(len(snps))
print "MNPs = % d" %(len(mnps))

# Insertion and deletion size distribution
insertion_sizes = [len(x) for x in insertions]
deletions_sizes = [len(x) for x in deletions]
r('insertions <- %s' % IntVector(insertion_sizes).r_repr())
r('deletions <- %s' % IntVector(deletions_sizes).r_repr())


r('''
pdf(file = "plots/indel_size_distri.pdf", width = 3.15, height = 1.75, pointsize = 9)
par(mar = c(5,4.5,2,1)+0.1, cex = 0.8)
hist(insertions, xlab="Insertion size", breaks=100)
hist(deletions, xlab="Deletion size", breaks=100)
dev.off()
''')

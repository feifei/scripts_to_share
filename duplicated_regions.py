''' To characterize duplicated regions 
    From filtered Blast match
    GC distribution, size distribution, gene contents, 
    copy number (MIN len)
'''

import argparse
import rpy2.robjects as robjects
from rpy2.robjects import *
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC, GC123
import re
import numpy

def get_seq_dict(scf_file):
    ''' Return sequence dictionary and total chromosome length '''
    
    seq_d = defaultdict(Seq)
    tot_chr_len = 0
    with open(scf_file, 'r') as scf_fh:
        for record in SeqIO.parse(scf_fh, 'fasta'):
            seq_d[record.id] = record.seq
            if record.id.startswith("chr"):
                tot_chr_len += len(record.seq)
    return seq_d, tot_chr_len


def filter_string(s):
    ''' Convert Artemis specific string to user friendly string'''
    
    return s.replace("%2C", ',').replace("+", " ").replace("%28", "(").replace("%09", "").\
             replace("%29", ")").replace("%3D", "=").replace("%2B", "+")


def get_annotation_dict(gff_file):
    ''' Return genes dictionary '''
    
    genes_d = defaultdict(list)
    with open(gff_file, 'r') as gff_fh:
        for line in gff_fh:
            arr = line.strip().split("\t")
            scfid = arr[0]
            feature = arr[2]
            if feature != "gene":
                continue
            start, end = map(int, arr[3:5])
            strand = arr[6]
            attrib = arr[8]
            geneid = re.match("ID=(.*?);", attrib).groups()[0]
            description = filter_string(re.search("label=.*?_(.*?);", attrib).groups()[0])
            genes_d[scfid].append([start, end, strand, geneid, description])
    return genes_d


def merge_intervals(intervals):
    ''' Merge overlapping intervals ''' 
    
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged
    

def get_gc_len(seq_d):
    ''' Return list of GC, sizes, and sequences of given sequences'''
    
    gcs = []
    lens = []
    seqs = []
    for scfid, l in seq_d.iteritems():
        for seq in l:
            gcs.append(GC(seq))
            lens.append(len(seq))
            seqs.append(seq)
    return gcs, lens, seqs


def get_gene_gc_len(genes_d):
    ''' Return list of gene GC and GC3, gene sizes, and detailed gene informations'''
    gcs = []
    gc3s = []
    lens = []
    info = []
    for scfid, l in genes_d.iteritems():
        for geneid, gene_start, gene_end, strand, description, gene_seq in l:
            gc = GC(gene_seq)
            gc3 = GC123(gene_seq)[3]
            length = len(gene_seq)
            gcs.append(gc)
            gc3s.append(gc3)
            lens.append(length)
            info.append([scfid, geneid, strand, length, gc, gc3, gene_start, gene_end, strand, description])
    return gcs, gc3s, lens, info


parser = argparse.ArgumentParser(description='Chracteraize duplicated regions')
parser.add_argument('blastn_file')
parser.add_argument('scf_fasta')
parser.add_argument('gff_annotation')
parser.add_argument('outfile')
parser.add_argument('--len_cutoff', dest='len_cutoff', default=1000, type=int, help='Length cutoff for the matches')
parser.add_argument('--ident_cutoff', dest='ident_cutoff', default=95, type=int, help='Identity cutoff for the matches')


#outfile = "../analysis/duplicated_regions/dupliciated_regions.txt"

args = parser.parse_args()
blastn_file = args.blastn_file
scf_file = args.scf_fasta
gff_file = args.gff_annotation
outfile = args.outfile
len_cutoff = args.len_cutoff
ident_cutoff = args.ident_cutoff
print(len_cutoff, ident_cutoff)



seq_d, tot_chr_len = get_seq_dict(scf_file)
genes_d = get_annotation_dict(gff_file)

pos_per_scf = defaultdict(list)
nr_within_chr = 0
nr_across_chr = 0
with open(blastn_file, 'r') as blastn_fh:
    for line in blastn_fh:
        line = line.strip()
        arr = line.split()
        scf1, scf2 = arr[0:2]
        idy = float(arr[2])
        aln_len = int(arr[3])
        s1, e1, s2, e2 = map(int, arr[6:10])
        if s1 > e1:
            s1, e1 = e1, s1
        if s2 > e2:
            s2, e2 = e2, s2
        
        if scf1.startswith('chr') and scf2.startswith('chr') and idy >= ident_cutoff and aln_len >= len_cutoff:
            if scf1 == scf2 and s1 == s2 and e1 and e2:
                # match to itself
                continue
            
            pos_per_scf[scf1].append([s1, e1])
            pos_per_scf[scf2].append([s2, e2])
            
            if scf1 == scf2:
                nr_within_chr += 1
            else:
                nr_across_chr += 1
    



# Getting duplicated regions/genes and non-duplicated regions
dup_seq_d = defaultdict(list)
non_dup_seq_d = defaultdict(list)
dup_genes_d = defaultdict(list)
for scfid, arr in pos_per_scf.iteritems():
    non_overlap_arr = merge_intervals(arr)
    seq = seq_d[scfid]
    genes = genes_d[scfid]
    len_seq = len(seq)
    non_dup_start = 1
    print(scfid)
    for s, e in non_overlap_arr:
        non_dup_seq = seq[non_dup_start - 1 : s]
        non_dup_seq_d[scfid].append(non_dup_seq)            
        dup_seq = seq[s - 1 : e]
        dup_seq_d[scfid].append(dup_seq)
        non_dup_start = e + 1
        
        # Check the gene contents in the duplicated regions
        genes_sorted_by_lower_bound = sorted(genes, key=lambda tup: tup[0])
        for gene_start, gene_end, strand, geneid, description in genes_sorted_by_lower_bound:
            if gene_start >= s and gene_end <= e:
                gene_seq = seq[gene_start-1:gene_end]
                dup_genes_d[scfid].append([geneid, gene_start, gene_end, strand, description, gene_seq])        
    
    non_dup_seq = seq[non_dup_start - 1 : len_seq]
    non_dup_seq_d[scfid].append(non_dup_seq)


dup_gc_l, dup_len_l, dup_seq_l = get_gc_len(dup_seq_d)
non_dup_gc_l, non_dup_len_l, non_dup_seq_l = get_gc_len(non_dup_seq_d)
dup_gene_gc_l, dup_gene_gc3_l, dup_gene_len_l, dup_info = get_gene_gc_len(dup_genes_d)


# Print out the genes in duplicated regions in a file
with open(outfile, 'w') as outh:
    for arr in dup_info:
        outh.write("\t".join(map(str, arr)) + "\n")


non_dup_genes_d = defaultdict(list)
for scfid, l in genes_d.iteritems():
    for gene_start, gene_end, strand, geneid, description in l:
        flag = False
        for dup_geneid, dup_gene_start, dup_gene_end, dup_strand, dup_description, dup_gene_seq in dup_genes_d[scfid]:
            if geneid == dup_geneid:
                flag = True
                break
        
        if not flag:
            non_dup_genes_d[scfid].append([geneid, gene_start, gene_end, strand, description, seq_d[scfid][gene_start-1:gene_end]])

non_dup_gene_gc_l, non_dup_gene_gc3_l, non_dup_gene_len_l, _ = get_gene_gc_len(non_dup_genes_d)


tot_dup_len = sum(dup_len_l)
tot_non_dup_len = sum(non_dup_len_l)
sum_len = tot_dup_len + tot_non_dup_len

print("Total duplicated region length: %d " %tot_dup_len)
print("Maximum length of duplicated region: %d" %(max(dup_len_l)))
print("Minimum length of duplicated region: %d" %(min(dup_len_l)))
print("Median length of duplicated region: %d" %(numpy.median(dup_len_l)))

print("Total non-duplicated region length: %d " %tot_non_dup_len)
print("Sum of duplicated + non-duplicated region length: %d " %sum_len)
print("Total genome chr length: %d" %tot_chr_len)
print("Percent of duplicated regions: %.2f%% " %(tot_dup_len/float(tot_chr_len)*100) )
print("%d within chromosome matches; %d across chromosome matches" %(nr_within_chr, nr_across_chr))

mean_GC = GC("".join(map(str, seq_d.values())))
print("Whole genome GC %.2f; Avg GC of duplicated genes %.2f; Avg GC of non-duplicated genes %.2f; \
Avg GC3 of duplicated genes %.2f; Avg GC3 of non-duplicated genes %.2f " \
       %(mean_GC, numpy.mean(dup_gene_gc_l), numpy.mean(non_dup_gene_gc_l), \
         numpy.mean(dup_gene_gc3_l), numpy.mean(non_dup_gene_gc3_l)))
print("Median GC of duplicated genes %.2f; Median GC of non-duplicated genes %.2f; \
Median GC3 of duplicated genes %.2f; Median GC3 of non-duplicated genes %.2f " \
       %(numpy.median(dup_gene_gc_l), numpy.median(non_dup_gene_gc_l), \
         numpy.median(dup_gene_gc3_l), numpy.median(non_dup_gene_gc3_l)))

mean_dup_GC = GC("".join(map(str, dup_seq_l)))
mean_non_dup_GC = GC("".join(map(str, non_dup_seq_l)))
print("Avg GC of duplicated regions %.2f; Avg GC of non-duplicated regions %.2f " \
       %(mean_dup_GC, mean_non_dup_GC))
       

non_dup_gc_l = FloatVector(non_dup_gc_l).r_repr()
non_dup_len_l = IntVector(non_dup_len_l).r_repr()
dup_gc_l = FloatVector(dup_gc_l).r_repr()
dup_len_l = IntVector(dup_len_l).r_repr()
mean_GC = IntVector([mean_GC]).r_repr()
non_dup_gene_gc_l = FloatVector(non_dup_gene_gc_l).r_repr()
non_dup_gene_gc3_l = FloatVector(non_dup_gene_gc3_l).r_repr()
non_dup_gene_len_l = IntVector(non_dup_gene_len_l).r_repr()
dup_gene_gc_l = FloatVector(dup_gene_gc_l).r_repr()
dup_gene_gc3_l = FloatVector(dup_gene_gc3_l).r_repr()
dup_gene_len_l = IntVector(dup_gene_len_l).r_repr()


r('non_dup_gc_l <- %s' % (non_dup_gc_l))
r('non_dup_len_l <- %s' % (non_dup_len_l))
r('dup_gc_l <- %s' % (dup_gc_l))
r('dup_len_l <- %s' % (dup_len_l))
r('mean_GC <- %s' % (mean_GC))
r('non_dup_gene_gc_l <- %s' % (non_dup_gene_gc_l))
r('non_dup_gene_gc3_l <- %s' % (non_dup_gene_gc3_l))
r('non_dup_gene_len_l <- %s' % (non_dup_gene_len_l))
r('dup_gene_gc_l <- %s' % (dup_gene_gc_l))
r('dup_gene_gc3_l <- %s' % (dup_gene_gc3_l))
r('dup_gene_len_l <- %s' % (dup_gene_len_l))


r('''
library(ggplot2)
pdf(file = "plots/duplicated_regions.pdf" )
non_dup.df <- data.frame(length=non_dup_len_l, GC=non_dup_gc_l)
dup.df <- data.frame(length=dup_len_l, GC=dup_gc_l)
non_dup.df$type <- "non_dup"
dup.df$type <- "dup"
merged <- rbind(non_dup.df, dup.df)
g1 <- ggplot(subset(merged, length < 20000), aes(length, fill = type)) + geom_density(alpha = 0.2) 
g2 <- ggplot(merged, aes(GC, fill = type)) + geom_density(alpha = 0.2)
g3 <- ggplot(subset(merged, length < 20000), aes(length, fill = type)) + geom_histogram(alpha = 0.5, bins = 100, aes(y = ..density..), position = 'identity')
g4 <- ggplot(merged, aes(GC, fill = type)) + geom_histogram(alpha = 0.5, bins = 100, aes(y = ..density..), position = 'identity')
g2 <- g2 + geom_vline(xintercept=mean_GC) + ggtitle("GC% of non-duplicated region vs. duplicated regions")
g4 <- g4 + geom_vline(xintercept=mean_GC) + ggtitle("GC% of non-duplicated region vs. duplicated regions")
print(g1)
print(g2)
print(g3)
print(g4)

non_dup_gene.df <- data.frame(length=non_dup_gene_len_l, GC=non_dup_gene_gc_l, GC3 = non_dup_gene_gc3_l)
dup_gene.df <- data.frame(length=dup_gene_len_l, GC=dup_gene_gc_l, GC3 = dup_gene_gc3_l)
non_dup_gene.df$type <- "non_dup"
dup_gene.df$type <- "dup"
merged_gene <- rbind(non_dup_gene.df, dup_gene.df)
g1 <- ggplot(merged_gene, aes(length, fill = type)) + geom_density(alpha = 0.2)
g2 <- ggplot(merged_gene, aes(GC, fill = type)) + geom_density(alpha = 0.2)
g3 <- ggplot(merged_gene, aes(length, fill = type)) + geom_histogram(alpha = 0.5, bins = 100, aes(y = ..density..), position = 'identity')
g4 <- ggplot(merged_gene, aes(GC, fill = type)) + geom_histogram(alpha = 0.5, bins = 100, aes(y = ..density..), position = 'identity')
g2 <- g2 + geom_vline(xintercept=mean_GC) + ggtitle("GC% of non-duplicated genes vs. duplicated genes")
g4 <- g4 + geom_vline(xintercept=mean_GC) + ggtitle("GC% of non-duplicated genes vs. duplicated genes")
g5 <- ggplot(merged_gene, aes(GC3, fill = type)) + geom_histogram(alpha = 0.5, bins = 100, aes(y = ..density..), position = 'identity')
print(g1)
print(g2)
print(g3)
print(g4)
print(g5)

dev.off()
''')


    
        
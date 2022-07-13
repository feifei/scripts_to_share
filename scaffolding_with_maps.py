'''
    Merge contigs into scaffolds based on optical maps
    based on placement report from MapSolver 
    (Manually sorted single mapping)
    Output unmapped contigs in a separate fasta file
'''

from optparse import OptionParser
from collections import defaultdict
import os, re
from Bio import SeqIO
from Bio.Seq import Seq
from lcs import *
from Bio.Emboss.Applications import WaterCommandline

def water_align(query_seq, target_seq):
    water_cline = WaterCommandline(asequence="asis:" + query_seq,
                                   bsequence="asis:" + target_seq,
                                   aformat="simple",
                                   gapopen=10,
                                   gapextend=0.5,
                                   outfile='stdout'
                                   )
    out_data, err = water_cline()
    return out_data



parser = OptionParser("usage: %prog contigs placement.reports")

(options, args) = parser.parse_args()

if len(args) != 2:
    parser.error("incorrect number of arguments")

contigs_file = args[0]
placement_file = args[1]

basename, extension = os.path.splitext(contigs_file)

scaffold_file = basename + ".scf.fasta"

contigs_d = defaultdict(Seq)
with open(contigs_file, 'r') as inh:
    for record in SeqIO.parse(inh, 'fasta'):
        ctg_id = record.id.split("|")[0]
        contigs_d[ctg_id] = record.seq


scf_d = defaultdict(Seq)
prev_real_chr_end = 0
prev_chr = ""
with open(placement_file, 'r') as ph:
    for line in ph:
        if not re.match("Spironucleus", line):
            continue
        
        arr = line.split("\t")
        if len(arr) != 7 :
            continue
        
        chr_id, chr_start, chr_end, ctg_id, ctg_start, ctg_end, ori = arr
        chr_start, chr_end, ctg_start, ctg_end = int(chr_start), int(chr_end), int(ctg_start), int(ctg_end)
        ori = int(ori)
        ctg_id = ctg_id.split("|")[1]
        m = re.search("Chromosome (\d)", chr_id)
        chr_id = m.groups()[0] if m else 1
        print chr_id, ctg_id, ori
        ctg_len = len(contigs_d[ctg_id])
        
        
        # Extend the chr start, end to the whole length scaffold
        if ori == 1:
            real_chr_start = chr_start - (ctg_start - 1)
            real_chr_end = chr_end + (ctg_len - ctg_end)
        else:
            real_chr_start = chr_start - (ctg_len - ctg_end)
            real_chr_end = chr_end + (ctg_start - 1)
        
                
        
        
        if chr_id != prev_chr: 
            # New chromosome
            # Gap to the chromosome end:
            gap_size = real_chr_start - 1

            if gap_size > 0:
                Ns = Seq("N" * gap_size)
                scf_d[chr_id] = Ns
            else:
                scf_d[chr_id] = Seq("")
                print "Sequence longer than the mapper", chr_id, ctg_id

            scf_d[chr_id] +=  contigs_d[ctg_id] if ori == 1 else contigs_d[ctg_id].reverse_complement()
                        
        elif real_chr_start >= prev_real_chr_end:
            # Real gap size
            gap_size = real_chr_start - prev_real_chr_end - 1
            Ns = Seq("N" * gap_size)
            scf_d[chr_id] += Ns
            scf_d[chr_id] +=  contigs_d[ctg_id] if ori == 1 else contigs_d[ctg_id].reverse_complement()
            
        elif real_chr_start < prev_real_chr_end:
            # Two matching scaffolds overlap
            # Extract the overlapped part, and align with water
            # if ctg_id == prev_ctg:
            #     # Speciell case, where the contig is splited with the first part reverse complemented
            #
            #     continue
            overlap_size = prev_real_chr_end - real_chr_start + 1
            overlap1_seq = scf_d[chr_id][-overlap_size :]
            overlap2_seq = contigs_d[ctg_id][0: overlap_size + 1] if ori == 1 else contigs_d[ctg_id].reverse_complement()[0: overlap_size + 1]
            aln_out = water_align(overlap1_seq, overlap2_seq)
            
            out_split = aln_out.split("\n")
            identity = float(re.search("\((.*)\)", out_split[25]).group(1).replace("%", ""))
            seq1_start = int(re.search("asis\s+(\d+)", out_split[32]).group(1))
            seq2_start = int(re.search("asis\s+(\d+)", out_split[34]).group(1))
            seq1_end = int(re.search("asis\s+(\d+)", out_split[-8]).group(1))
            seq2_end = int(re.search("asis\s+(\d+)", out_split[-6]).group(1))
            print identity, overlap_size, seq1_start, seq1_end, seq2_start, seq2_end
            # if identity < 85:
            #     print aln_out
            scf_d[chr_id] = scf_d[chr_id][0: -(overlap_size - seq1_end)] + contigs_d[ctg_id][seq2_end:] if ori == 1 else \
                            scf_d[chr_id][0: -(overlap_size - seq1_end)] + contigs_d[ctg_id].reverse_complement()[seq2_end:]
            
        else:
            print "Else, finns det?"
            
        
        prev_real_chr_end = real_chr_end
        prev_chr = chr_id
        prev_ctg = ctg_id



with open(scaffold_file, 'w') as outh:
    for chr_id, seq in sorted(scf_d.items()):
        print >>outh, ">chr%s" %chr_id
        print >>outh, seq
        
        

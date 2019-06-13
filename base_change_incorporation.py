'''
    Incorporation of changes from pre-selected BWA mapping pileup results
    Indel and SNPs. 
    Modify the fasta file, and the corresponding gff file
'''
import os
import argparse
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


parser = argparse.ArgumentParser(description='Update gff file from pre-selected pileup base changes results, \
                                              and optioanlly changes in the reference file too')
parser.add_argument('base_changes_manual')
parser.add_argument('reference_file')
parser.add_argument('gff_file')
parser.add_argument('--fix_reference', action='store_true')


args = parser.parse_args()
base_changes_file = args.base_changes_manual
reference_file = args.reference_file
gff_file = args.gff_file
new_reference_file = reference_file + ".changes"
new_gff_file = gff_file + ".changes"
fix_reference = args.fix_reference # default False


global_pos_d = defaultdict(int)
ref_seq_d = defaultdict(str)
tot_len = 0
ordered_scfids = []
with open(reference_file, 'r') as inh:
    for record in SeqIO.parse(inh, 'fasta'):
        scfid = record.id
        ref_seq_d[scfid] = record.seq
        global_pos_d[scfid] = tot_len
        tot_len += len(record.seq)
        ordered_scfids.append(scfid)



def insert_str(string, insert, index):
    # Insert after the given index
    return string[:index] + insert + string[index:]


def delete_str(string, delete, index):
    return string[:(index - len(delete) + 1)] + string[index + 1:]


def replace_str(string, old, new, index):
    return string[:(index - len(old))] + new + string[index:]
    

new_ref_seq_d = defaultdict(str)
changes_d = defaultdict(list)
extra = 0
prev_scfid = ""
snp_count = 0
indel_count = 0
with open(base_changes_file, 'r') as inh:
    for line in inh:
        scfid, orig_pos, orig_base, new_base = line.strip().split()[0:4]
        orig_pos = int(orig_pos)
        global_pos = global_pos_d[scfid] + orig_pos
        
        if prev_scfid != scfid:
            extra = 0
            if prev_scfid:
                ref_seq_d[prev_scfid] = scf_seq
            scf_seq = ref_seq_d[scfid]
        
            
        if len(orig_base) > 3 or len(new_base) > 3:
            print scfid, pos, global_pos, orig_base, new_base
            continue
        
        idx_pos = int(orig_pos) + extra    
        if orig_base == ".":
            # Insertion
            scf_seq = insert_str(scf_seq, new_base, idx_pos)
            print "Insertion", scfid, global_pos, new_base, scf_seq[idx_pos - 5 : idx_pos + 5]
            extra += len(new_base)
            change_type = "insertion"
            indel_count += 1
        elif new_base == ".":
            # Deletion
            scf_seq = delete_str(scf_seq, orig_base, idx_pos)
            print "Deletion", scfid, global_pos, orig_base, scf_seq[idx_pos - 5 : idx_pos + 5]
            extra -= len(orig_base)
            change_type = "deletion"
            indel_count += 1
        else:
            # SNPs
            scf_seq = replace_str(scf_seq, orig_base, new_base, idx_pos)
            print "SNP", scfid, global_pos, orig_base, new_base, scf_seq[idx_pos - 5 : idx_pos + 5]
            extra = extra - len(orig_base) + len(new_base)
            change_type = "SNP"
            snp_count += 1
        
        changes_d[scfid].append([orig_pos, extra, change_type])
        prev_scfid = scfid


ref_seq_d[scfid] = scf_seq
records = []
for scfid in ordered_scfids:
    records.append(SeqRecord(ref_seq_d[scfid], id=scfid, description = ""))


if fix_reference:
    # Write out the new sequences in the fasta format
    with open(new_reference_file, 'w') as outh:
        SeqIO.write(records, outh, 'fasta')


def calculate_new_pos(pos_info, input_pos):
    # Given extra pos info to modify the start, end in the gff annotation file
    new_pos = 0
    prev_extra = 0
    if input_pos > pos_info[-1][0]:
        new_pos = input_pos + pos_info[-1][1]
    else:
        for pos, extra, _ in pos_info:
            if input_pos < pos:
                new_pos = input_pos + prev_extra
                continue
            prev_extra = extra
    return new_pos

    

def is_in_gene(pos_info, start, end):
    # To check if the gene has been modified based on changes
    in_gene_count = 0
    change_types = []
    for change_pos, _, change_type in pos_info:
        if start <= change_pos and end > change_pos:
            # Gene affected
            in_gene_count += 1
            change_types.append(change_type)
            
    return [in_gene_count, change_types]


snp_gene_count = 0
indel_gene_count = 0
with open(gff_file, 'r') as gff_inh, open(new_gff_file, 'w') as gff_outh:
    for line in gff_inh:
        line = line.strip()
        if len(line) == 0:
            continue
        arr = line.split("\t")
        scfid = arr[0]
        start, end = map(int, arr[3:5])
        if scfid in changes_d:
            new_start = calculate_new_pos(changes_d[scfid], start)
            new_end = calculate_new_pos(changes_d[scfid], end)
        
            arr[3] = new_start
            arr[4] = new_end
            
            in_gene_count, change_types = is_in_gene(changes_d[scfid], start, end)
            if in_gene_count:
                for change_type in change_types:                    
                    if change_type == "SNP":
                        snp_gene_count += 1
                        print "SNP",
                    else:
                        indel_gene_count += 1
                        print "Indel",
                print "\n" + line
                
        
        print >>gff_outh, "\t".join(map(str, arr))


print "Total snp changes %d, and %d in genes; Total indel changes %d, and %d in genes" %(snp_count, snp_gene_count, indel_count, indel_gene_count)
import sys
import subprocess
from Bio import SeqIO


def trainGlimmerHMM(mfasta_file, exon_file, args):
    ''' train GlimmerHMM with core genes
        core genes are prepared as a tab delimited exon files.
        trained results are written into a train directory
    '''
    command = ['trainGlimmerHMM', mfasta_file, exon_file] \
              + sum(map(list, zip(args.keys(), args.values())), [])
    cline = subprocess.Popen(command, stdout = subprocess.PIPE)
    cline.stdout.close()
    assert(cline.wait() == 0)


def glimmerhmm(genome_file, train_dir_for_genome, args):
    ''' Run glimmerHMM to predict genes/exons '''
    
    command = ['glimmerhmm', genome_file, train_dir_for_genome] \
            + sum(map(list, zip(args.keys(), args.values())), [])
    cline = subprocess.Popen(command, stdout = subprocess.PIPE)
    cline.stdout.close()
    assert(cline.wait() == 0)

home_dir = "glimmerhmm/"
genome_file = "ssk.cns.fa"
exon_file = "training_for_glimmerhmm.cds"
train_dir = home_dir + "training_dir4"
conf_file = train_dir + "/train_0_100.cfg"
gff_output = home_dir + "glimmerhmm.prediction4.gff"


trainGlimmerHMM(genome_file, exon_file, {'-d': train_dir,'-n':'150', '-v':'50'})

# Modify the train config file to add onlytga=1
with open(conf_file, 'a') as conf_h:
    print >>conf_h, "BoostSgl 10"
    print >>conf_h, "onlytga 1"

# Since glimmerhmm takes only one sequence at a time,
# Loop through the scaffolds, and predict separately, then concatenate 
with open(gff_output, 'w') as gff_h:
    print >>gff_h, "##gff-version 3"
    
    with open(genome_file, 'r') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seqid = record.id
            seqfile = seqid + ".fa"
            with open(seqfile, 'w') as seqh:
                SeqIO.write(record, seqh, 'fasta')
    
            tmp_output = home_dir + seqid + ".gff"
            glimmerhmm(seqfile, train_dir, {'-o': tmp_output, '-g':'', '-v':''})
            
            with open(tmp_output, 'r') as tmp_h:
                for i, line in enumerate(tmp_h):
                    if i != 0:
                        gff_h.write(line)
            
            subprocess.check_call('rm %s %s' % (seqfile, tmp_output), shell = True)


    
    

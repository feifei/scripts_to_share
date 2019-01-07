import sys
import subprocess

def prodigal(in_file, out_file, args):
    cmd = ['prodigal',
            '-i', in_file, '-o', out_file] \
            + sum(map(list, zip(args.keys(), args.values())), [])
    cline = subprocess.Popen(cmd, stdout = subprocess.PIPE)
    cline.stdout.close()
    assert(cline.wait() == 0)
    return out_file

genome_file = "ssk.cns.fa"

'''
So providing training file doesn't really work in prodigal
Or impossible, since the training file is a binary file
It does train on the provided genome though
training_file = "data/core_genes.cns.fa"
Used training_file as input, and write out the training file
which is then used to predict the genome file
'''

program = "prodigal"
outfile = program + "/" + program + ".prediction.gff"
prodigal(genome_file, outfile, {'-g':'6', '-f':"gff"})

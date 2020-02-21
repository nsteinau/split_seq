#!/usr/bin/env python

import sys
import argparse

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-F', '--inputFastqF', required=True, help='Input Fastq File')
args = parser.parse_args()

   
# Define eprint function to print to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

id_used=[]

#trimmded_filename="/media/myxfs/splitseq/SPLiT-Seq_demultiplexing/combined_run_results/AGCATTCGAACCGAGA-TGGAACAAATCCA-ACAAGCTAGTGGCC_tx_trimmed.fastq"
with open(args.inputFastqF, "r") as infile:
#with open(trimmded_filename,"r") as infile: 
    line_ct = 0
    for line in infile:
        if (line_ct % 4 == 0):
            line=line.split("_")
            id_used.append(line[0])
        line_ct += 1

file = args.inputFastqF.split("/")
#file = trimmded_filename.split("/")
file = file[-1]
file = file[:-17]
file_fastq = str("/media/xfs/splitseq/SPLiT-Seq_demultiplexing/combined_results_run4/umi_results/" + file + ".fastq")
file_umi = str("/media/xfs/splitseq/SPLiT-Seq_demultiplexing/combined_results_run4/" + file + ".umi")


umi=[]
# Read in input .fastq and store lines as a python dictionary
with open(file_fastq, "r") as infile2:
    line_ct = 0
    for line in infile2:
        if (line_ct % 4 == 0):
            line=line.split("_")
            if line[0] in id_used:
                umi.append(line[5])
        line_ct += 1

with open(file_umi, "w+") as f:
	for key in umi:
		f.writelines(key)

sys.exit()

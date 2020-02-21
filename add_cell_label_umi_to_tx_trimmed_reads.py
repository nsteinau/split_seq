#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 16:49:48 2020

@author: jinsong
"""
import argparse
from itertools import islice
parser = argparse.ArgumentParser()
parser.add_argument('-F', '--inputFastqF', required=True, help='Input Fastq File')
args = parser.parse_args()

   
# Define eprint function to print to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

            
umi_file=args.inputFastqF
#umi_file="/media/myxfs/splitseq/SPLiT-Seq_demultiplexing/combined_results_run2/AGCATTCGTTCACGCA-TTCACGCAATCCA-GTCTGTCAGTGGCC.umi"
cell_id=umi_file.split("/")
cell_id=cell_id[-1]
cell_id=cell_id[:-4]
cell_id
cell_id=cell_id.replace("-","")
umi=open(umi_file,"r")
umi_id=[]
for line in umi:
    line=line.strip()
    umi_id.append(line)
cdna_file=umi_file[:-4] + "_tx_trimmed.fastq"
idx=0
processed_reads=[]
with open(cdna_file,"r") as f:
    while True:
        read = list(islice(f, 4))
        if not read:
            break
        id = read[0].splitlines()[0]
        read[0]=id+"_"+cell_id + "_" + umi_id[idx]+"\n"
        for item in read:
            processed_reads.append(item)
        idx = idx + 1
process_file=umi_file[:-49] + "tx_trimmed_for_star/" + umi_file[-49:] + ".fastq"
with open(process_file,"w+") as wf:
    wf.writelines(processed_reads)

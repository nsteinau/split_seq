#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import sys
import splitseq_utilities
from datetime import datetime
from itertools import islice
import numpy as np
import argparse


# In[ ]:

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--inputFastq', required=True, help='Input Fastq File')
parser.add_argument('-o', '--outDir', required=True, help='Directory to store demultiplexed fastq files')
args = parser.parse_args()

# Define eprint function to print to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# In[ ]:


os.chdir("/media/xfs/splitseq/SPLiT-Seq_demultiplexing/")
read_stat = {}
statistics = {}
buffers = {}
round1barcodeDictionary = {}
round2barcodeDictionary = {}
round3barcodeDictionary = {}
round1_odt_hex={}
possibleChars = ['A', 'C', 'G', 'T', 'N',""]
LENGTH = "LENGTH"
HYPHEN = "-"
FASTQ_EXTENSION = ".fastq"


# In[ ]:


class FastQRead(object):
    def __init__(self):
        self.round1barcodeMatches = set()
        self.round2barcodeMatches = set()
        self.round3barcodeMatches = set()
        self.data = []

class BarcodeMatch(object):
    def __init__(self, barcode, startIndex, endIndex,rank):
        self.barcode = barcode
        self.startIndex = startIndex
        self.endIndex = endIndex
        self.rank=rank
    def __eq__(self, other):
        if isinstance(other, BarcodeMatch):
            return self.barcode == other.barcode
        return False
    def __hash__(self):
        return hash(self.barcode)

class ReadIndex(object):
    def __init__(self, index):
        self.index = index

def barcodesToDictionary(barcodeDictionary, barcodeFile, errors):

    file = open(barcodeFile, "r")

    for barcode in file:
        barcode = barcode.strip()
        if LENGTH not in barcodeDictionary:
            barcodeDictionary[LENGTH] = []
            barcodeDictionary[LENGTH].append(len(barcode))
            for i in range(-errors, errors + 1):
                if i != 0:
                    barcodeDictionary[LENGTH].append(len(barcode) - i)
        splitseq_utilities.addToDictionarySet(barcodeDictionary, barcode, barcode)
        if errors != 0:
            barcodeVariantsToDictionary(barcodeDictionary, barcode, barcode, 0, 0, errors, set())
    file.close()


def barcodeVariantsToDictionary(barcodeDictionary, barcode, variant, index, depth, errors, memo):

    #Termination Clause
    if depth >= errors:
        return

    #DP Termination Clause
    memoItem = variant + "-" + str(depth)
    if memoItem in memo:
        return

    memo.add(memoItem)

    for i in range(index, len(variant)):
        currentChar = variant[i]
        for possibleChar in possibleChars:
            if possibleChar != "":
                newVariant = variant[:i] + possibleChar + variant[i:]
                splitseq_utilities.addToDictionarySet(barcodeDictionary, newVariant, barcode)
                barcodeVariantsToDictionary(barcodeDictionary, barcode, newVariant, i + 1, depth + 1, errors, memo)
            if currentChar != possibleChar:
                newVariant = variant[:i] + possibleChar + variant[i+1:]
                splitseq_utilities.addToDictionarySet(barcodeDictionary, newVariant, barcode)
                barcodeVariantsToDictionary(barcodeDictionary, barcode, newVariant, i + 1, depth + 1, errors, memo)


def crawlFastQ(fastqr, outputdir, targetMemory, granularity):

    startTime = datetime.now()
    counter = 0

    with open(fastqr) as f:

        while True:
            read = FastQRead()
            read.data = list(islice(f, 4))
            if not read.data:
                break
            line = read.data[1]
            readIndex = ReadIndex(len(line) - 1)
            read.round1barcodeMatches = crawlSegments(line, readIndex, round1barcodeDictionary, round1barcodeDictionary[LENGTH])
            if read.round1barcodeMatches:

                read.round2barcodeMatches = crawlSegments(line, readIndex, round2barcodeDictionary, round2barcodeDictionary[LENGTH])
                if read.round2barcodeMatches:

                    read.round3barcodeMatches = crawlSegments(line, readIndex, round3barcodeDictionary, round3barcodeDictionary[LENGTH])

            addToBuffers(read)
            real = FastQRead()
            counter += 1
            if (counter % granularity) == 0:
                print("Analyzed [" + "{:,}".format(counter) + "] reads in [{}]".format(datetime.now() - startTime))
                memoryUsage = splitseq_utilities.analyzeMemoryUsage()
                if memoryUsage > targetMemory:
                    print("\tCurrent memory [{}] exceeded target memory [{}]. Flushing buffers...".format(splitseq_utilities.bytesToDisplay(memoryUsage), splitseq_utilities.bytesToDisplay(targetMemory)))
                    splitseq_utilities.flushBuffers(outputdir, buffers)
    print("Analyzed [" + "{:,}".format(counter) + "] reads in [" + str(datetime.now() - startTime) + "]")
    splitseq_utilities.flushBuffers(outputdir, buffers)


def crawlSegments(line, readIndex, barcodeDictionary, barcodeLengths):

    results = set()
    nextIndex = 0
    for barcodeLength in barcodeLengths:

        index = readIndex.index
        while index > 10 + barcodeLength - 1:

            segment = line[index-barcodeLength:index]
            if segment in barcodeDictionary:

                if segment in barcodeDictionary[segment]:

                    if not results:
                        results.add(BarcodeMatch(segment, index-barcodeLength, index,0))
                        if nextIndex < index - barcodeLength + 1:
                            nextIndex = index - barcodeLength + 1
                    else:
                        for item in results:

                            if item.rank > 0:
                                item.barcode = segment
                                item.startIndex = index-barcodeLength
                                item.endIndex = index
                                item.rank = 0
                                if nextIndex < index - barcodeLength + 1:
                                    nextIndex = index - barcodeLength + 1

                else:
                    if not results:

                        if len(barcodeDictionary[segment]) == 1:
                            for barcode in barcodeDictionary[segment]:
                                results.add(BarcodeMatch(barcode, index-barcodeLength, index,1))
                                if nextIndex < index - barcodeLength + 1:
                                    nextIndex = index - barcodeLength + 1

            index = index - 1


    readIndex.index = nextIndex
    return results

def addToBuffers(read):
    read_id = read.data[0].splitlines()[0]
    if read.round1barcodeMatches and read.round2barcodeMatches and read.round3barcodeMatches:
        for round1barcodeMatch in read.round1barcodeMatches:
            for round2barcodeMatch in read.round2barcodeMatches:
                for round3barcodeMatch in read.round3barcodeMatches:
                    # Filter out false positives that didn't find barcodes in the proper order
                    if round2barcodeMatch.endIndex > round1barcodeMatch.startIndex:
                        continue
                    if round3barcodeMatch.endIndex > round2barcodeMatch.startIndex:
                        continue
                    if round1barcodeMatch.barcode in round1_odt_hex.keys():
                        round1barcodeMatch.barcode = round1_odt_hex[round1barcodeMatch.barcode]
                    bufferKey = round1barcodeMatch.barcode + HYPHEN + round2barcodeMatch.barcode + HYPHEN + round3barcodeMatch.barcode + FASTQ_EXTENSION
                    read.data[0] = read.data[0].splitlines()[0] + "_" + round1barcodeMatch.barcode + "_" + str(round1barcodeMatch.startIndex) + "_" + str(round1barcodeMatch.endIndex) + "\n"                   
                    for line in read.data:
                        splitseq_utilities.addToDictionaryList(buffers, bufferKey, line)

                    #Track total number of blocks added to each file
                    if bufferKey in statistics:
                        statistics[bufferKey] += 1
                    else:
                        statistics[bufferKey] = 1
                    if read_id in read_stat:
                        read_stat[read_id] += 1
                    else:
                        read_stat[read_id] = 1


# In[ ]:


minreads=10
round1barcodes="Round1_barcodes_new5.txt"
round2barcodes="Round2_barcodes_new4.txt"
round3barcodes="Round3_barcodes_new4.txt"
#fastqr="kasumi_split_seq/R_2019_12_20_11_43_08_user_S5-00516-12-SPLIT_Seq_Kasumi_1_mixed_SPLITSeq.fastq"
#fastqr="/media/xfs/splitseq/SPLiT-Seq_demultiplexing/kasumi_split_seq/combined_split_seq.fastq"
fastqr=args.inputFastq
#fastqr="/media/xfs/splitseq/SPLiT-Seq_demultiplexing/kasumi_split_seq/combined_test.fastq"
errors=1
#outputdir="combined_results_run3"
outputdir=args.outDir
#outputdir="combined_test_run2"
targetMemory=40000000000
granularity=1000000


# In[ ]:


round1_odt_hex={}
odt_hex=[]
with open(round1barcodes,"r") as f:
    for line in f:
        line=line.splitlines()
        odt_hex.append(line[0])
for i in range(48):
    round1_odt_hex[odt_hex[i+48]]=odt_hex[i]


# In[ ]:


preprocessingStartTime = datetime.now()
barcodesToDictionary(round1barcodeDictionary, round1barcodes, errors)
barcodesToDictionary(round2barcodeDictionary, round2barcodes, errors)
barcodesToDictionary(round3barcodeDictionary, round3barcodes, errors)
print("Pre-processing data structures completed in [{}]".format(datetime.now() - preprocessingStartTime))
print("\tCurrent memory [{}]...".format(splitseq_utilities.bytesToDisplay(splitseq_utilities.analyzeMemoryUsage())))
splitseq_utilities.createDirectory(outputdir)
splitseq_utilities.clearFilesMatchingFilter(outputdir, lambda f: True)
crawlFastQ(fastqr, outputdir, targetMemory, granularity)

print("# of results files [{}]".format(splitseq_utilities.countFilesMatchingFilter(outputdir, lambda f: f.endswith("fastq"))))
splitseq_utilities.clearFilesMatchingFilter(outputdir, lambda f: f.endswith("fastq") and statistics[f] < minreads)
print("# of results files [{}]".format(splitseq_utilities.countFilesMatchingFilter(outputdir, lambda f: f.endswith("fastq"))))
two_cells = []
one_cell = []
for i in read_stat.keys():
    if read_stat[i]==2:
        two_cells.append(i)
    else:
        one_cell.append(i)
print("number of reads assigned to at least two cells = ", two_cells)

from itertools import islice
import argparse
parser = argparse.ArgumentParser()
def hammingtondistance(str1,str2):
    l=len(str1)
    hd=0
    for x in range(l):
        if str1[x:x+1]!=str2[x:x+1]:
            hd=hd+1
    return hd
parser.add_argument("--file", "-f", type=str, required=True)
args = parser.parse_args()
adaptor="CTGTCTCTTATACACATCT"
class FastQRead(object):
    def __init__(self):
        self.round1barcodeMatches = set()
        self.round2barcodeMatches = set()
        self.round3barcodeMatches = set()
        self.data = []
with open(args.file) as f:
    processed_reads=[]
    while True:
        read=FastQRead()
        read.data = list(islice(f, 4))
        if not read.data:
            break
        line=read.data[1]
        line=line.splitlines()[0]
        idx=len(str(line))
        while idx > 19:
            d=hammingtondistance(line[(idx-19):idx],adaptor)
            if d > 2:
                idx=idx-1
                continue
            else:
                if line[idx-19] == adaptor[0]:
                    read.data[1]=read.data[1][:(idx-19)]  + "\n"
                    read.data[3]=read.data[3][:(idx-19)]  + "\n"
                    break
                else:
                    if line[idx-18] == adaptor[1]:
                        read.data[1]=read.data[1][:(idx-18)]  + "\n"
                        read.data[3]=read.data[3][:(idx-18)]  + "\n"
                        break
                    else:
                        read.data[1]=read.data[1][:(idx-17)]  + "\n"
                        read.data[3]=read.data[3][:(idx-17)]  + "\n"
                        break
        if len(read.data[1]) > 20:
            processed_reads.append(read.data)
if len(processed_reads) > 0:
    with open(args.file[:-6] + "_" + "trimmed.fastq","w") as tx_file:
        for x in processed_reads:
            tx_file.writelines(x)

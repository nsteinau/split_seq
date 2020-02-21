from itertools import islice
import argparse
parser = argparse.ArgumentParser()                                               

parser.add_argument("--file", "-f", type=str, required=True)
args = parser.parse_args()
read_tx=[]
read=[]
with open(args.file) as f:
        while True:
            read = list(islice(f, 4))
            if not read:
                break
            id = read[0]
            id=id.splitlines()[0]
            id=id.split("_")
            id=int(id[-1])
            read[1]=read[1][id:]
            read[3]=read[3][id:]
            if read[1] =="\n":
                read[1] = "N\n"
                read[3] = "N\n"
            read_tx.append(read)
with open(args.file[:-6] + "_" + "tx.fastq","w") as tx_file:
    for x in read_tx:
        tx_file.writelines(x)

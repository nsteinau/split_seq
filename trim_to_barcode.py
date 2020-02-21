import pandas as pd
from snakemake.shell import shell

shell("seqkit locate -j {threads} -m 1 -f {input.bc3} | cut -f1,2,3,4,5  > barcode3hits.txt")
shell("seqkit locate -j {threads} -m 1 -f {input.bc1} | cut -f1,6  > barcode1hits.txt")

BC3 = pd.read_fwf("barcode3hits.txt")

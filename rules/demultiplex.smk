rule make_bc_fasta:
    input:
        bc1 = config["bc1"],
        bc2 = config["bc2"],
        bc3 = config["bc3"]
    output:
        bc1 = temp("barcodes/bc1.fa"),
        bc2 = temp("barcodes/bc2.fa"),
        bc3 = temp("barcodes/bc3.fa")
    shell:
        """
        cat {input.bc1} | awk '{{printf ">%s\\n%s\\n",$0,$0}}' > {output.bc1}
        cat {input.bc2} | awk '{{printf ">%s\\n%s\\n",$0,$0}}' > {output.bc2}
        cat {input.bc3} | awk '{{printf ">%s\\n%s\\n",$0,$0}}' > {output.bc3}
        """

rule find_bcs:
    input:
        reads = config["raw_reads"],
        bc1 = "barcodes/bc1.fa",
        bc2 = "barcodes/bc2.fa",
        bc3 = "barcodes/bc3.fa"
    output:
        bc1 = "barcodes/bc1hits.bed",
        bc3 = "barcodes/bc3hits.bed"
    threads: 24
    shell:
        """
        seqkit locate -j {threads} -m 1 --bed -M -P -f {input.bc3} {input.reads} | cut -f1,2 > {output.bc3}
        seqkit locate -j {threads} -m 1 --bed -M -P -f {input.bc1} {input.reads} | cut -f1,3  > {output.bc1}
        """

rule extract_bcs:
    input:
        bc1 = "barcodes/bc1hits.txt",
        bc3 = "barcodes/bc3hits.txt"
    output:
        barcodes = "barcodes/positions.txt"
    run:
        BC3 = pd.read_csv(input.bc3,sep="\t")
        BC1 = pd.read_csv(input.bc1,sep="\t")
        pos = pd.merge(BC3,BC1,how='inner',on='seqID')
        pos['bc_length'] = pos.end - pos.start
        pos.sort_values('bc_length',ascending=False).drop_duplicates(subset='seqID')

#! /bin/bash

# Step 1: demultiplex the single FASTQ file into many FASTQ files belonging to single cells + collapse oligoDT and random-Hexmer barcodes

python split_seq_demultiplex_ver0.23_collapse_oligodt_randomhex.24_collapse_oligoDT_randomhex.py -f $1 -o $2

# Step 2: process the reads to remove the 5' end that contains barcode1/2/3.

#extarct_cdna.py

export OUTPUT_DIR=$1

ls $OUTPUT_DIR | grep \.fastq$ | parallel -j 30 -k "python extarct_cdna.py -f $OUTPUT_DIR/{}"


# Step 3: Trim 3' portion of the reads to remove Tn5 univeral adaptor, also remove reads shorter than 20 bps

#trim_split_seq.py


ls $OUTPUT_DIR | grep _tx\.fastq$ | parallel -j 30 -k "python trim_split_seq.py -f $OUTPUT_DIR/{}"


#Step 4: extract umi located upstream of barcode 1 and downstream of Illumina TrueSeq and put it at the end of seq ID.
 
#umi_tools

mkdir $OUTPUT_DIR/umi_results

ls $OUTPUT_DIR | grep \.fastq$ | grep -v tx | parallel -j 30 -k "umi_tools extract --stdin=$OUTPUT_DIR/{} --bc-pattern=CAGACGTGTGCTCTTCCGATCTNNNNNNNNNN --log=processed.log --stdout=$OUTPUT_DIR/umi_results/{}"

#Step 5: extract Umi into a separte files for each read contained in the output file of trimmed cDNA reads from Step 3

#align_kallisto_kasumi_1_split_seq_trimmed_test.py

ls $OUTPUT_DIR | grep \tx_trimmed.fastq$ | parallel -j 30 -k "python align_kallisto_kasumi_1_split_seq_trimmed_test.py -F $OUTPUT_DIR/{}"

#Step 6: Add cell labela and umi into the seq id for every read in fastq files outputted from Step 3 and Step 5

#add_cell_label_umi_to_tx_trimmed_reads_v2_add_import_sys.py

mkdir $OUTPUT_DIR/tx_trimmed_for_star

ls $OUTPUT_DIR | grep \.umi$ | parallel -j 30 -k "python add_cell_label_umi_to_tx_trimmed_reads.py -F $OUTPUT_DIR/{}"


#Step 7: aggregate all fastq files into one fastq file

cd $OUTPUT_DIR/tx_trimmed_for_star

mkdir all_cdna

ls | grep \.umi\.fastq | parallel -j 30 -k "cat {}" > all_cdna/all_cdna_tx_trimmed_all.fastq

STAR --runThreadN 10 --readFilesIn all_cdna_tx_trimmed_all.fastq  --outFileNamePrefix ./all_tx_cdna_mm2 --genomeDir /media/data_analysis5/apps/STAR/index/GRCh38.p13.genome  --alignIntronMax 20000  --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2

fc=/media/data_analysis5/apps/subread-2.0.0-Linux-x86_64/bin/featureCounts

$fc -F GTF -a /media/data_analysis5/apps/STAR/index/GRCh38.p13.genome/gencode.v32.annotation.gtf  -o gene_assigned -R BAM all_tx_cdna_mm2Aligned.sortedByCoord.out.bam  -T 20  -M

samtools sort all_tx_cdna_mm2Aligned.sortedByCoord.out.bam.featureCounts.bam  -o all_tx_cdna_mismatch2_assigned_sorted.bam

samtools index all_tx_cdna_mismatch2_assigned_sorted.bam

conda activate bioinformatics
umi_tools count --wide-format-cell-counts --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I all_tx_cdna_mismatch2_assigned_sorted.bam -S counts.tsv.gz

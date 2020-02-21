import pandas as pd
from snakemake.utils import validate


##### load config and sample sheets #####

configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")


##### target rules #####

rule all:
    input:"barcodes/bc1hits.bed","barcodes/bc3hits.bed"



##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
#singularity: "docker://continuumio/miniconda3"






##### load rules #####
include: "rules/demultiplex.smk"

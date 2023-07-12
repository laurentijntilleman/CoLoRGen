#!/bin/bash

# parameters

## overall variables
# constants
number_treats=40

# paths to the different files
nanopore_sequencing_files='fastq_pass'  #raw fastq files
analyse_folder='./' #working directory
ref_folder="ref/" # reference folder

# paths to the reference genomes
reference_genome='GRCh38.fa'
minimap2_index='GRCh38.fa.mmi'
reference_gtf='GRCh38.gtf'


# edit paths to external software tools
export PATH=samtools-1.11:$PATH
export PATH=bcftools-1.11:$PATH
export PATH=htslib-1.11:$PATH

hap_py='hap.py' #path to the hap.py script

# gene information
#CYP2D6: GRCh38: NC_000022.11:c42130810-42126499
#CYP2D7: GRCh38: NC_000022.11:c42144483-42139576
# targeted region
gene_json='cyp2d6_hybrids.json'
star_alleles="${ref_folder}CYP2D6.NC_000022.11.haplotypes.tsv"

target_start=42121165
target_stop=42149371
target_chr='chr22'

## variant calling parameters
# medaka model
model_m='r941_min_hac_g507'

log_level='debug'

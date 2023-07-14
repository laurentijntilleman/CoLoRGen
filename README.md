CoLoRGen
========

CoLoRGen: comprehensive long read genotyping pipeline.
This tool generates consensus sequence of the targeted regio per allele
and returns vcf files for each gene that is found. If star-alleles are
present, the correct star-allles are assigned to the genes.

This is CoLoRGen v2 in nextflow


Installation
------------

CoLoRGen can be run by downloading these git repository. CoLoRGen uses conda to install all the packages and nextflow to run the pipeline.

Usage
-----

Running CoLoRGen is easy. With the command below the full pipeline is run.

    nextflow run colorgen.nf

The following paramters need to be adjust in the `nextflow.config` file:

    index          = 'GRCh38.fa'   # fasta file of the reference genome
    reference_gtf  = 'GRCh38.gtf'  # gtf file with the exon regions
    reads          = 'fastq_pass'  # raw fastq files
    threads        = 40            # number of threads
    input_bed      = 'target.bed'  # the targeted region
    gene_json      = 'cyp2d6_hybrids.json' # json file with gene info and hybrid info
    ref_folder     = 'CYP2D6.NC_000022.11.haplotypes.tsv' # star allele nomenclature downloaded from PharmVar
    output_dir     = 'output'      # path to the output folder
    model_m        = 'r941_min_hac_g507' # model for medaka consensus
    clair3_model   = '/bin/models/r941_prom_sup_g5014' # model for clair3 variant calling 


In the git repository, there is a example of the json file with the gene info


Output
------

The output is structured in different folder. All the intermediate files are
kept for debugging en visualisation proposes.

## Output structure

    .
    ├── ref                     # folder with genome and star-allele references
    ├── haplotype
    │   ├── consensus_1         # first consensus based on the structural variants
    │   ├── consensus_1_breakpoints.txt # break points defining the structural variants
    │   ├── consensus_2         # second consensus constructed with the reads that mapped on the alleles defined by medaka
    │   ├── consensus_2_H1      # consensus 2 folder for H1
    │   ├── consensus_2_H2      # consensus 2 folder for H2
    │   ├── consensus_3         # third consensus constructed with all the reads
    │   ├── consensus_3_genes   # target genes mapped on consensus sequence
    │   ├── consensus_3_H1      # consensus 3 folder for H1
    │   ├── consensus_3_H2      # consensus 3 folder for H2
    │   ├── coverage.txt        # depth on the target sequence
    │   ├── haplotypes_H1_gene_*_*.vcf # vcf file of the genes on allele 1
    │   ├── haplotypes_H2_gene_*_*.vcf # vcf file of the genes on allele 2
    │   ├── haplotypes_haplotypes.txt  # overview of the found star-alleles
    │   ├── haplotypes.vcf      # vcf file of all the found variants
    │   ├── medaka_round_0      # medaka output of the input data
    │   ├── select_reads_taged  # reads split by haplotypes defined by medaka
    ├── log
    │   ├── error.log           # error log
    │   └── stout.log           # standard log
    ├── mapped                  # mapped data
    └── raw_reads               # raw reads

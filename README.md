CoLoRGen
========

CoLoRGen: comprehensive long read genotyping pipeline.
This tool generates consensus sequence of the targeted regio per allele
and returns vcf files for each gene that is found. If star-alleles are
present, the correct star-allles are assigned to the genes.


Installation
------------

CoLoRGen can be installed by downloading these git repository and run the
install script.

Using this method requires the user to provide several binaries:

 * [samtools](https://github.com/samtools/samtools),
 * [minimap2](https://github.com/lh3/minimap2),
 * [tabix](https://github.com/samtools/htslib), and
 * [bgzip](https://github.com/samtools/htslib)

`samtools/bgzip/tabix` version 1.11 and `minimap2` version 2.18 are recommended
as these are those used in development of CoLoRGen.

Most of the script use python3, therefore python3 need to be installed on the
system. Also the virtualenv package need to be present as CoLoRGen runs
everything in a virtual environment.

Usage
-----

Running CoLoRGen is easy. With the command below the full pipeline is run.

    CoLoRGen -p <parameters.sh>

The following paramters need to be adjust in the `parameters.sh` file:

    number_treats=40
    nanopore_sequencing_files='fastq_pass'  # raw fastq files
    analyse_folder='./' # working directory
    ref_folder="ref/" # reference folder
    reference_genome='GRCh38.fa' # fasta file of the reference genome
    minimap2_index='GRCh38.fa.mmi' # minimap2 index file
    reference_gtf='GRCh38.gtf' # gtf file with the exon regions
    # edit paths to external software tools
    export PATH=samtools-1.11:$PATH
    export PATH=bcftools-1.11:$PATH
    export PATH=htslib-1.11:$PATH
    hap_py='hap.py' #path to the hap.py script
    gene_json='cyp2d6_hybrids.json' # json file with gene info and hybrid info
    star_alleles="${ref_folder}CYP2D6.NC_000022.11.haplotypes.tsv" # star allele nomenclature downloaded from PharmVar
    # position of the targeted region
    target_start=42121165
    target_stop=42149371
    target_chr='chr22'
    # medaka model
    model_m='r941_min_hac_g507' # model for medaka variant calling and consensus
    log_level='debug' # log level

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

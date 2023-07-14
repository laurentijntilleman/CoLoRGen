// ~/tools/nextflow/nextflow run colorgen.nf -resume --reads /data/projects/20220705_AS_pharmacogenes/20230405_NA12878_analysis_sup_select/pass_select.fastq.gz --index /data/ensembl_genomes/Hsa/release107/Homo_sapiens.GRCh38.dna.primary_assembly.fa -with-conda true


// index reference
process index_fasta {
  conda 'samtools'
  
  input:
    path fasta
  output:
    tuple path("$fasta"), path("${fasta}.fai")
  
  """
  samtools faidx $fasta
  """
 }

 process build_gene{
  conda 'colorgen_scripts.yml'

  input:
    path gene_json
    tuple path(index), path(index_fai)

  output:
    path "output/*"

  """
    mkdir output
    build_gene.py -i $gene_json -o output -r $index
  """


 }


// align the reads to the reference geneome
process align {
  conda 'minimap2 samtools'

  input:
    tuple path(index), path(index_fai)
    path reads
  output:
    tuple path('sorted.bam'), path('sorted.bam.bai')
    
  """
  minimap2 -t $params.threads -a $index $reads -2 | samtools view -@ $params.threads -Sb | samtools sort -@ $params.threads > sorted.bam
  samtools index -@ $params.threads sorted.bam
  """
}


// Create a channel of BED intervals
bed_intervals = Channel.fromPath(params.input_bed)
  .splitCsv(header:false, sep: '\t')
  .map { row -> [ row[0], row[1] as int, row[2] as int , row[3]] }


// Define a process to split the BAM file
process split_bam {
  conda 'samtools'
  
  input:
    tuple path(bam_file), path(x)
    tuple val(chrom), val(start), val(end), val(region) 
  output:
    tuple path("${region}.bam"), path("${region}.bam.bai"), val(region)
  
  """
  samtools view -@ $params.threads -b $bam_file ${chrom}:${start}-${end} > "${region}.bam"
  samtools index -@ $params.threads "${region}.bam"

  """
}

process phase_reads {
  conda 'clair3.yml'
  
  input:
    tuple path(region_bam), path(region_bam_bai), val(region)
    tuple path(index), path(index_fai)
    
  output:
    tuple path("${region}_taged.bam"), path("${region}_taged.bam.bai"), val(region)
  
  """
  run_clair3.sh --bam_fn=$region_bam  --ref_fn=$index --threads=$params.threads --platform="ont" --model_path=\${CONDA_PREFIX}${params.clair3_model} --output=$region
  whatshap phase -o ${region}_phased.vcf --reference=$index "${region}/merge_output.vcf.gz" $region_bam --distrust-genotypes --ignore-read-groups --indels
  bgzip ${region}_phased.vcf
  tabix ${region}_phased.vcf.gz
  whatshap haplotag --ignore-read-groups -o ${region}_taged.bam --reference ${index} ${region}_phased.vcf.gz $region_bam
  samtools index ${region}_taged.bam
  """
}

process split_reads {
  conda 'colorgen_scripts.yml'
  
  input:
    tuple path(region_taged_bam), path(region_taged_bam_bai), val(region)
    tuple val(chrom), val(start), val(end), val(region) 
    
  output:
    tuple path("${region}_taged_split.fastq.gz"), path("${region}_taged_H1.fastq.gz"), path("${region}_taged_H2.fastq.gz"), path("${region}_taged.fastq.gz")
  
  """
  get_second_aligments.py -i $region_taged_bam -o ${region}_taged -s $start -t $end
  cat ${region}_taged_H*.fastq | gzip > ${region}_taged.fastq.gz
  gzip ${region}_taged_split_sec.fastq
  gzip ${region}_taged_split.fastq
  gzip ${region}_taged_H*.fastq
  """
}

process align_sec {
  conda 'minimap2 samtools'

  input:
    tuple path(index), path(index_fai)
    tuple path(reads_taged_split_fastq), path(reads_taged_H1_fastq), path(reads_taged_H2_fastq), path(reads_taged_fastq_gz)
    tuple val(chrom), val(start), val(end), val(region) 

  output:
    tuple path("${region}_sorted.bam"), path("${region}_sorted.bam.bai")
    
  """
  minimap2 -t $params.threads -a $index $reads_taged_split_fastq -2 | samtools view -@ $params.threads -Sb | samtools sort -@ $params.threads > ${region}_sorted.bam
  samtools index -@ $params.threads ${region}_sorted.bam
  """
}

process build_consensus {
  conda 'colorgen_scripts.yml'

  publishDir "${params.output_dir}", pattern: "*.fasta"

  input:
    tuple path(sorted_bam), path(sorted_bam_bai)
    tuple path(index), path(index_fai)
    tuple val(chrom), val(start), val(end), val(region) 
  output:
    tuple path("${region}_consensus_1_breakpoints.txt"), path("${region}_consensus_1_H1.fasta"), path("${region}_consensus_1_H2.fasta")

  """
  touch "${region}_${chrom}.txt"
  build_consensus_sequence.py -i $sorted_bam -o ${region}_consensus_1 -r $index -c $chrom -s ${start} -t ${end}
  
  """

}

process polish_consensus_H1 {
  conda 'colorgen_scripts_medaka.yml'

  input:
    tuple path(consensus_1_breakpoints), path(consensus_1_H1_fasta), path(consensus_1_H2_fasta)
    tuple path(taged_split_fastq), path(taged_H1_fastq), path(taged_H2_fastq), path(taged_fastq)
    val model_m
    tuple val(chrom), val(start), val(end), val(region) 

  output:
    tuple path("${region}_consensus_H1.bam"), path("${region}_consensus_H1.bam.bai"), path("${region}_consensus_2_H1/consensus.fasta")

  """
  minimap2 -t $params.threads -a $consensus_1_H1_fasta $taged_H1_fastq -2 | samtools view -Sb | samtools sort > ${region}_consensus_H1.bam
  samtools index ${region}_consensus_H1.bam
  medaka_consensus -i $taged_H1_fastq -d ${region}_consensus_1_H1.fasta -o ${region}_consensus_2_H1 -t $params.threads -m ${model_m}
  """
}

process polish_consensus_H2 {
  conda 'colorgen_scripts_medaka.yml'

  input:
    tuple path(consensus_1_breakpoints), path(consensus_1_H1_fasta), path(consensus_1_H2_fasta)
    tuple path(taged_split_fastq), path(taged_H1_fastq), path(taged_H2_fastq), path(taged_fastq)
    val model_m
    tuple val(chrom), val(start), val(end), val(region) 

  output:
    tuple path("${region}_consensus_H2.bam"), path("${region}_consensus_H2.bam.bai"), path("${region}_consensus_2_H2/consensus.fasta")

  """
  minimap2 -t $params.threads -a $consensus_1_H2_fasta $taged_H2_fastq -2 | samtools view -Sb | samtools sort > ${region}_consensus_H2.bam
  samtools index ${region}_consensus_H2.bam
  medaka_consensus -i $taged_H2_fastq -d ${region}_consensus_1_H2.fasta -o ${region}_consensus_2_H2 -t $params.threads -m ${model_m}
  """
}

process align_concensus_2 {
  conda 'colorgen_scripts.yml'

  input:
    tuple path(consensus_H1_bam), path(consensus_H1_bam_bai), path("${region}_consensus_2_H1/consensus.fasta")
    tuple path(consensus_H2_bam), path(consensus_H2_bam_bai), path("${region}_consensus_2_H2/consensus.fasta")
    tuple path(taged_split_fastq), path(taged_H1_fastq), path(taged_H2_fastq), path(taged_fastq)
    tuple val(chrom), val(start), val(end), val(region) 

  output:
    tuple path("${region}_C2_clean.bam"), path("${region}_C2_clean.bam.bai"), path("${region}_C2_clean_H1.fastq"), path("${region}_C2_clean_H2.fastq")
  
  """
  cat ${region}_consensus_2_H1/consensus.fasta ${region}_consensus_2_H2/consensus.fasta > ${region}_C2.fasta

  minimap2 -t $params.threads -a ${region}_C2.fasta $taged_fastq -2 | samtools view -Sb | samtools sort > ${region}_sorted.bam
  samtools index ${region}_sorted.bam
  clean_bam.py -i ${region}_sorted.bam -o ${region}_C2_clean
  samtools index ${region}_C2_clean.bam
  """

}

process polish_consensus2_H1 {
  conda 'colorgen_scripts_medaka.yml'

  input:
    tuple path(clean_bam), path(clean_bam_bai), path(clean_H1_fastq), path(clean_H2_fastq)
    tuple path(consensus_bam), path(consensus_bam_bai), path(consensus_fasta)
    val model_m
    tuple val(chrom), val(start), val(end), val(region) 

  output:
    path "${region}_C3_H1/consensus.fasta"

  """
  medaka_consensus -i $clean_H1_fastq -d $consensus_fasta -o ${region}_C3_H1 -t $params.threads -m ${model_m}
  """
}

process polish_consensus2_H2 {
  conda 'colorgen_scripts_medaka.yml'

  input:
    tuple path(clean_bam), path(clean_bam_bai), path(clean_H1_fastq), path(clean_H2_fastq)
    tuple path(consensus_bam), path(consensus_bam_bai), path(consensus_fasta)
    val model_m
    tuple val(chrom), val(start), val(end), val(region) 

  output:
    path "${region}_C3_H2/consensus.fasta"

  """
  medaka_consensus -i $clean_H2_fastq -d $consensus_fasta -o ${region}_C3_H2 -t $params.threads -m ${model_m}
  """
}

process align_concensus_3 {
  conda 'colorgen_scripts.yml'

  publishDir "${params.output_dir}", pattern: "${region}*"

  input:
    tuple path(clean_bam), path(clean_bam_bai), path(clean_2_H1_fastq), path(clean_2_H2_fastq)
    path "${region}_C3_H1/consensus.fasta"
    path "${region}_C3_H2/consensus.fasta"
    tuple val(chrom), val(start), val(end), val(region) 
    
  output:
    tuple path("${region}_C3_clean.bam"), path("${region}_C3_clean.bam.bai"), path("${region}_consensus.fasta")
 
 """
  cat ${region}_C3_H1/consensus.fasta ${region}_C3_H2/consensus.fasta > ${region}_consensus.fasta
  cat $clean_2_H1_fastq $clean_2_H2_fastq | gzip > clean_2.fastq.gz

  minimap2 -t $params.threads -a ${region}_consensus.fasta clean_2.fastq.gz -2 | samtools view -Sb | samtools sort > sorted.bam
  samtools index sorted.bam
  clean_bam.py -i sorted.bam -o ${region}_C3_clean
  samtools index ${region}_C3_clean.bam
  """

}

process align_genes {
  conda 'bowtie2'

  publishDir "${params.output_dir}", pattern: "${region}*"

  input:
    tuple path(clean_bam), path(clean_bam_bai), path(fasta)
    tuple val(chrom), val(start), val(end), val(region)
    path gene_folder_fasta
  
  output:
    tuple path("${region}_sorted.bam"), path("${region}_sorted.bam.bai")

  """
  bowtie2-build $fasta sample
  bowtie2 -a -x sample -p $params.threads -f genes.fasta | samtools view -Sb | samtools sort > ${region}_sorted.bam
  samtools index ${region}_sorted.bam

  """
}

process calculate_probs{
  conda 'colorgen_scripts_medaka.yml'

  input:
    tuple path(clean_bam), path(clean_bam_bai), path(fasta)
    val model_m
    tuple val(chrom), val(start), val(end), val(region)

  output:
    path "${region}_probs.hdf"


  """
    medaka consensus $clean_bam ${region}_probs.hdf --model $model_m --batch_size 100  --threads $params.threads
  """
}

process annotate_variants {
  conda 'colorgen_scripts_medaka.yml'

  publishDir 'output'

  input:
    tuple path(clean_bam), path(clean_bam_bai), path(fasta)
    tuple path(gene_sorted_bam), path(gene_sorted_bam_bai)
    path probs
    tuple val(chrom), val(start), val(end), val(region) 
    path star_alleles
    path gene_json

  output:
    path "${region}_*"

  """
    variant_caller.py $fasta $gene_sorted_bam $probs $star_alleles $gene_json ${region} $params.reference_gtf
  """
}

workflow {
  index             = index_fasta(params.index)
  targeted_genes    = build_gene(params.gene_json,index)
  bam_sorted        = align(index,Channel.fromPath(params.reads))
  region_bam        = split_bam(bam_sorted,bed_intervals)
  phase_bam         = phase_reads(region_bam,index)
  reads_split       = split_reads(phase_bam,bed_intervals)
  bam_sec           = align_sec(index,reads_split,bed_intervals)
  consensus_1       = build_consensus(bam_sec,index,bed_intervals)
  consensus_2_H1    = polish_consensus_H1(consensus_1,reads_split,params.model_m,bed_intervals)
  consensus_2_H2    = polish_consensus_H2(consensus_1,reads_split,params.model_m,bed_intervals)
  consensus_2_clean = align_concensus_2(consensus_2_H1,consensus_2_H2,reads_split,bed_intervals)
  consensus_3_H1    = polish_consensus2_H1(consensus_2_clean,consensus_2_H1,params.model_m,bed_intervals)
  consensus_3_H2    = polish_consensus2_H2(consensus_2_clean,consensus_2_H2,params.model_m,bed_intervals)
  consensus_3       = align_concensus_3(consensus_2_clean,consensus_3_H1,consensus_3_H2,bed_intervals)
  probs             = calculate_probs(consensus_3,params.model_m,bed_intervals)
  aligned_genes     = align_genes(consensus_3,bed_intervals,targeted_genes)
  annotate_variants(consensus_3,aligned_genes,probs,bed_intervals,params.star_alleles,params.gene_json)
}


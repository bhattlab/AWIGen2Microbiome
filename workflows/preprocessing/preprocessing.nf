#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
 * pipeline input parameters -- This should go into a config file before long
 */
params.projectdir = "/labs/asbhatt/dmaghini/tools/nextflow_tutorial/indata"
params.reads = "/labs/asbhatt/dmaghini/tools/nextflow_tutorial/indata/*{1,2}.fq.gz"


log.info """\
    METAGENOMIC PREPROCESSING - NF PIPELINE
    ===================================
    raw.reads        : ${params.reads}
    host.genome      : ${params.bwa_index_base}
    
    outdir reads     : ${params.outdir}
    """
    .stripIndent()

/* PREPROCESSING
 * Processes for fastqc, multiqc, and read processing, including
 * deduplication, trimming, human read removal
*/

process fastqc {
    publishDir params.outdir + "/stats/pre_fastqc", pattern: '*logs', mode: params.publish_mode
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "prefastqc_${sample_id}_logs", emit: logs
    path "counts_${sample_id}_raw.tsv", emit: stats

    script:
    """
    mkdir prefastqc_${sample_id}_logs
    fastqc -o prefastqc_${sample_id}_logs --extract -f fastq -q ${reads}
    readcount=\$(echo \$((\$(zcat ${reads[0]} | wc -l) / 2)))
    echo ${sample_id}"\tRaw\t"\$readcount > "counts_${sample_id}_raw.tsv"
    """
}

process multiqc {
    publishDir params.outdir + "/stats", mode: params.publish_mode

    input:
    path '*'

    output:
    path 'multiqc_pre_report.html'

    script:
    """
    multiqc --filename multiqc_pre_report.html .
    """
}

process deduplicate {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R*.fastq.gz"), emit: dedupreads
    path("counts_${sample_id}_dedup.tsv"), emit: dedupstats

    script:
    """
    hts_SuperDeduper -1 ${reads[0]} -2 ${reads[1]} -f ${sample_id} -F
    readcount=\$(echo \$((\$(zcat ${sample_id}_R1.fastq.gz | wc -l) / 2)))
    echo ${sample_id}"\tDeduplicate\t"\$readcount > "counts_${sample_id}_dedup.tsv"
    """
}


process trimgalore {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}*val_*.fq.gz"), emit: trimreads
    path("counts_${sample_id}_trim.tsv"), emit: trimstats

    script:
    """
    trim_galore --quality ${params.trimgalore_quality} --length ${params.trimgalore_min_read_length} \
        --paired ${reads[0]} ${reads[1]} --retain_unpaired

    zcat -f ${sample_id}_R1_unpaired_1.fq.gz ${sample_id}_R2_unpaired_2.fq.gz | pigz -b 32 > ${sample_id}_val_unpaired.fq.gz
    rm ${sample_id}_R1_unpaired_1.fq.gz ${sample_id}_R2_unpaired_2.fq.gz

    readcount_paired=\$(echo \$((\$(zcat ${sample_id}_R1_val_1.fq.gz | wc -l) / 2)))
    readcount_unpaired=\$(echo \$((\$(zcat ${sample_id}_val_unpaired.fq.gz | wc -l) / 4)))
    totalcount=\$(echo \$((\$readcount_paired + \$readcount_unpaired)))
    echo ${sample_id}"\tTrim\t"\$totalcount > "counts_${sample_id}_trim.tsv"
    """
}

process hostremoval {
    publishDir params.outdir  + "/preprocessed_reads", pattern: '*gz', mode: params.publish_mode

    input:
    tuple val(sample_id), path(reads)
    path host_genome_location
    val bwa_index_base

    output:
    tuple val(sample_id), path("${sample_id}_cleaned_*fastq.gz"), emit: hostremreads
    path("counts_${sample_id}_hostrem.tsv"), emit: hoststats

    script:
    """
    bwa mem ${host_genome_location}/${bwa_index_base} ${reads[0]} ${reads[1]} | \
        samtools fastq -t -T BX -f 4 -1 ${sample_id}_cleaned_1.fastq.gz -2 ${sample_id}_cleaned_2.fastq.gz -s ${sample_id}_cleanedtemp_singletons.fastq.gz -

        # run on unpaired reads
    bwa mem ${host_genome_location}/${bwa_index_base} ${reads[2]} | \
        samtools fastq -t -T BX -f 4  - > ${sample_id}_cleanedtemp_singletons2.fastq.gz

    # combine singletons
    zcat -f ${sample_id}_cleanedtemp_singletons.fastq.gz ${sample_id}_cleanedtemp_singletons2.fastq.gz | pigz > ${sample_id}_cleaned_orphans.fastq.gz
    rm ${sample_id}_cleanedtemp_singletons.fastq.gz ${sample_id}_cleanedtemp_singletons2.fastq.gz

    readcount_paired=\$(echo \$((\$(zcat ${sample_id}_cleaned_1.fastq.gz | wc -l) / 2)))
    readcount_unpaired=\$(echo \$((\$(zcat ${sample_id}_cleaned_orphans.fastq.gz | wc -l) / 4)))
    totalcount=\$(echo \$((\$readcount_paired + \$readcount_unpaired)))
    echo ${sample_id}"\tHostRemoved\t"\$totalcount > "counts_${sample_id}_hostrem.tsv"
    """

}

process postfastqc{
    publishDir params.outdir + "/stats/post_fastqc", mode: params.publish_mode
    input:
    tuple val(sample_id), path(reads)

    output:
    path "postfastqc_${sample_id}_logs"

    script:
    """
    mkdir postfastqc_${sample_id}_logs
    fastqc -o postfastqc_${sample_id}_logs --extract -f fastq -q ${reads}
    """
}

process postmultiqc {
    publishDir params.outdir + "/stats", mode: params.publish_mode

    input:
    path '*'

    output:
    path 'multiqc_postreport.html'

    script:
    """
    multiqc --filename multiqc_postreport.html .
    """
}

process aggregatereports {
  publishDir params.outdir + "/stats/read_counts", mode: params.publish_mode
  
  input:
  path rawstats
  path dedupstats
  path trimstats
  path hoststats

  output:
  path 'readcounts.tsv'

  script:
  """
  cat ${rawstats}  ${dedupstats} ${trimstats} ${hoststats} > readcounts.tsv
  """
}


/* ASSEMBLY
 * Runs megahit on paired and orphan reads, then uses QUAST to measure
 * assembly quality. Combines quast reports into a single report.
*/

process megahit {
  publishDir params.outdir + "/assembly", pattern: 'megahit_out/*contigs.fa', mode: params.publish_mode

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("megahit_out/${sample_id}.contigs.fa")

  shell:
  """
  megahit -1 ${reads[0]} -2 ${reads[1]} -r ${reads[2]} --out-prefix ${sample_id}
  """

}

process quast {
  input:
  tuple val(sample_id), path(contigs)

  output:
  tuple val(sample_id), path(quast)

  shell:
  """
  quast.py -o quast ${contigs} --fast
  """
}

workflow {
    /* FIXME
    The input should be a bit more customizable... will have to think about this in the future
    for example: what about other file endings? samples sequenced on several lanes?
    for the future, though
    */
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    // PREPROCESSING
    fastqc_ch = fastqc(read_pairs_ch)
    multiqc(fastqc_ch.logs.collect())
    deduplicated_ch = deduplicate(read_pairs_ch)
    trim_galore_ch = trimgalore(deduplicated_ch.dedupreads)
    host_remove_ch = hostremoval(trim_galore_ch.trimreads, params.host_genome_location, params.bwa_index_base)
    postfastqc_ch = postfastqc(host_remove_ch.hostremreads)
    postmultiqc(postfastqc_ch.collect())
    aggregatereports(fastqc_ch.stats.collect(), 
                     deduplicated_ch.dedupstats.collect(), 
                     trim_galore_ch.trimstats.collect(), 
                     host_remove_ch.hoststats.collect())
    
    // CLASSIFICATION
    
    // ASSEMBLY
    megahit_ch = megahit(host_remove_ch.hostremreads)
    //quast(megahit_ch)
    
    // BINNING

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Yay!\n" : "Oops .. something went wrong" )
}

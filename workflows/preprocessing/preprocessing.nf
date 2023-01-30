#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
 * pipeline input parameters -- This should go into a config file before long
 */
params.projectdir = "/labs/asbhatt/dmaghini/tools/nextflow_tutorial/indata"
params.reads = "/labs/asbhatt/dmaghini/tools/nextflow_tutorial/indata/*{1,2}.fq.gz"


log.info """\
    METAGENOMIC PREPROCESSING - NF PIPELINE
    =======================================
    raw.reads        : ${params.reads}
    host.genome      : ${params.bwa_index_base}
    
    outdir reads     : ${params.outdir_preprocessed_reads}
    outdir stats     : ${params.outdir_stats}
    """
    .stripIndent()

process fastqc {
    publishDir params.outdir_stats + '/pre_fastqc', mode: params.publish_mode
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
    publishDir params.outdir_stats, mode: params.publish_mode

    input:
    path '*'

    output:
    path 'multiqc_report_pre.html'

    script:
    """
    multiqc --filename multiqc_report_pre.html .
    """
}

process deduplicate {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R*.fastq.gz")

    /*
     * The functions gives a stats.log file which could be parsed to
     * get the number of emitted/duplicated reads (without counting the files)
     */
    script:
    """
    hts_SuperDeduper -1 ${reads[0]} -2 ${reads[1]} -f ${sample_id} -F
    """
}

/*
 * FIXME: still need to edit this to take in preprocessing values as params
 * and optionally to take in params for start/end trim in different positions
 */
process trimgalore {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}*val_*.fq.gz")

    script:
    """
    trim_galore --quality ${params.trimgalore_quality} --length ${params.trimgalore_min_read_length} \
      --paired ${reads[0]} ${reads[1]} --retain_unpaired

    zcat -f ${sample_id}_R1_unpaired_1.fq.gz ${sample_id}_R2_unpaired_2.fq.gz | pigz -b 32 > ${sample_id}_val_unpaired.fq.gz
    rm ${sample_id}_R1_unpaired_1.fq.gz ${sample_id}_R2_unpaired_2.fq.gz
    """
}

process hostremoval {
    publishDir params.outdir_preprocessed_reads, mode: params.publish_mode

    input:
    tuple val(sample_id), path(reads)
    path host_genome_location

    output:
    tuple val(sample_id), path("${sample_id}_cleaned_*fastq.gz")

    script:
    """
    bwa mem ${host_genome_location}/*.fa ${reads[0]} ${reads[1]} | \
        samtools fastq -t -T BX -f 4 -1 ${sample_id}_cleaned_1.fastq.gz -2 ${sample_id}_cleaned_2.fastq.gz -s ${sample_id}_cleanedtemp_singletons.fastq.gz -

        # run on unpaired reads
    bwa mem ${host_genome_location}/*.fa ${reads[2]} | \
        samtools fastq -t -T BX -f 4  - > ${sample_id}_cleanedtemp_singletons2.fastq.gz

    # combine singletons
    zcat -f ${sample_id}_cleanedtemp_singletons.fastq.gz ${sample_id}_cleanedtemp_singletons2.fastq.gz | pigz > ${sample_id}_cleaned_orphans.fastq.gz
    rm ${sample_id}_cleanedtemp_singletons.fastq.gz ${sample_id}_cleanedtemp_singletons2.fastq.gz
    """

}

process postfastqc{
    publishDir params.outdir_stats + '/post_fastqc', mode: params.publish_mode
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

/* FIXME still need to test this
*/
process postmultiqc {
    publishDir params.outdir_stats, mode: params.publish_mode

    input:
    path '*'

    output:
    path 'multiqc_report_post.html'

    script:
    """
    multiqc --filename multiqc_report_post.html .
    """
}


workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    fastqc_ch = fastqc(read_pairs_ch)
    multiqc(fastqc_ch.logs.collect())
    deduplicated_ch = deduplicate(read_pairs_ch)
    trim_galore_ch = trimgalore(deduplicated_ch)
    host_remove_ch = hostremoval(trim_galore_ch, params.host_genome_location)
    postfastqc_ch = postfastqc(host_remove_ch)
    postmultiqc(postfastqc_ch.collect())

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Yay!\n" : "Oops .. something went wrong" )
}


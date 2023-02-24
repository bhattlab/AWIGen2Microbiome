process fastqc {
    publishDir params.outdir + "/stats/pre_fastqc", mode: params.publish_mode
    tag "pre-FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "prefastqc_${sample_id}_logs/*"

    script:
    """
    mkdir prefastqc_${sample_id}_logs
    fastqc -o prefastqc_${sample_id}_logs -f fastq -q ${reads} -t ${task.cpus}
    """
}

process postfastqc{
    publishDir params.outdir + "/stats/post_fastqc", mode: params.publish_mode
    tag "post-FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "postfastqc_${sample_id}_logs"

    script:
    """
    mkdir postfastqc_${sample_id}_logs
    fastqc -o postfastqc_${sample_id}_logs -f fastq -q ${reads} -t ${task.cpus}
    """
}

process fastqc {
    publishDir params.outdir + "/stats/pre_fastqc/", mode: params.publish_mode, pattern: "*logs/*"
    tag "pre-FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "prefastqc_${sample_id}_logs/*", emit: prefastqc
    path "versions.yml", emit: versions

    script:
    """
    mkdir prefastqc_${sample_id}_logs
    fastqc -o prefastqc_${sample_id}_logs -f fastq -q ${reads} -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS

    """
}

process postfastqc{
    publishDir params.outdir + "/stats/post_fastqc", mode: params.publish_mode
    tag "post-FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "postfastqc_${sample_id}_logs", emit: postfastqc

    script:
    """
    mkdir postfastqc_${sample_id}_logs
    fastqc -o postfastqc_${sample_id}_logs -f fastq -q ${reads} -t ${task.cpus}
    """
}

/*
 * pipeline input parameters
 */
params.projectdir = "/labs/asbhatt/dmaghini/tools/nextflow_tutorial/indata"
params.reads = "/labs/asbhatt/dmaghini/tools/nextflow_tutorial/indata/*{1,2}.fq.gz"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"
log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */

process FASTQC {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    fastqc.sh "$sample_id" "$reads"
    """
}

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(fastqc_ch.collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}

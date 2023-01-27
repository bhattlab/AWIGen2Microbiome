/*
 * pipeline input parameters
 */
params.projectdir = "/labs/asbhatt/dmaghini/tools/nextflow_tutorial/indata"
params.reads = "/labs/asbhatt/dmaghini/tools/nextflow_tutorial/indata/*{1,2}.fq.gz"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"
params.bwa_index_base = "/labs/asbhatt/data/host_reference_genomes/hg19/hg19.fa"


log.info """\
    METAGENOMIC PREPROCESSING - N F   P I P E L I N E
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
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
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

process DEDUPLICATE {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R*.fastq.gz")

    script:
    """
    hts_SuperDeduper -1 ${reads[0]} -2 ${reads[1]} -p ${sample_id} -F -g
    """
}

/*
 * FIXME: still need to edit this to take in preprocessing values as params
 * and optionally to take in params for start/end trim in different positions
 */
process TRIMGALORE {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}*val_*.fq.gz")

    script:
    """
    trim_galore --quality 30 --length 60 --paired ${reads[0]} ${reads[1]} --retain_unpaired

    zcat -f ${sample_id}_R1_unpaired_1.fq.gz ${sample_id}_R2_unpaired_2.fq.gz | pigz -b 32 > ${sample_id}_val_unpaired.fq.gz
    rm ${sample_id}_R1_unpaired_1.fq.gz ${sample_id}_R2_unpaired_2.fq.gz
    """
}

process HOSTREMOVAL {
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_cleaned_*fastq.gz")

    script:
    """
    bwa mem /labs/asbhatt/data/host_reference_genomes/hg19/hg19.fa ${reads[0]} ${reads[1]} | \
        samtools fastq -t -T BX -f 4 -1 ${sample_id}_cleaned_1.fastq.gz -2 ${sample_id}_cleaned_2.fastq.gz -s ${sample_id}_cleanedtemp_singletons.fastq.gz -

        # run on unpaired reads
    bwa mem /labs/asbhatt/data/host_reference_genomes/hg19/hg19.fa ${reads[2]} | \
        samtools fastq -t -T BX -f 4  - > ${sample_id}_cleanedtemp_singletons2.fastq.gz

    # combine singletons
    zcat -f ${sample_id}_cleanedtemp_singletons.fastq.gz ${sample_id}_cleanedtemp_singletons2.fastq.gz | pigz > ${sample_id}_cleaned_orphans.fastq.gz
    rm ${sample_id}_cleanedtemp_singletons.fastq.gz ${sample_id}_cleanedtemp_singletons2.fastq.gz
    """

}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(fastqc_ch.collect())
    deduplicated_ch = DEDUPLICATE(read_pairs_ch)
    trim_galore_ch = TRIMGALORE(deduplicated_ch)
    trim_galore_ch.view()
    host_remove_ch = HOSTREMOVAL(trim_galore_ch)

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}

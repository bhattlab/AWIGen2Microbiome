/*
 * pipeline input parameters
 */
params.projectdir = "/labs/asbhatt/dmaghini/tools/nextflow_tutorial/indata"
params.reads = "/labs/asbhatt/dmaghini/tools/nextflow_tutorial/indata/*{1,2}.fq.gz"
params.multiqc = "$projectDir/multiqc"
params.outdir = "/labs/asbhatt/dmaghini/projects/awigen2/AWIGen2Microbiome/workflows/preprocessing/results"
params.bwa_index_base = "/labs/asbhatt/data/host_reference_genomes/hg19/hg19.fa"


log.info """\
    METAGENOMIC PREPROCESSING - N F   P I P E L I N E
    ===================================
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

process fastqc {
    publishDir params.outdir + "/stats", mode:'copy'
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "prefastqc_${sample_id}_logs", emit: logs
    path "counts_${sample_id}_raw.tsv", emit: stats

    script:
    """
    mkdir prefastqc_${sample_id}_logs
    fastqc -o prefastqc_${sample_id}_logs -f fastq -q ${reads}
    readcount=\$(echo \$((\$(zcat ${reads[0]} | wc -l) / 2)))
    echo ${sample_id}"\tRaw\t"\$readcount > "counts_${sample_id}_raw.tsv"
    """
}

process multiqc {
    publishDir params.outdir + "/stats", mode:'copy'

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
process trimgalore {
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

process hostremoval {
    publishDir params.outdir  + "/preprocessing", mode:'copy'

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

process postfastqc{
    publishDir params.outdir + "/stats"
    input:
    tuple val(sample_id), path(reads)

    output:
    path "postfastqc_${sample_id}_logs"

    script:
    """
    mkdir postfastqc_${sample_id}_logs
    fastqc -o postfastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

/* FIXME still need to test this
*/
process postmultiqc {
    publishDir params.outdir + "/stats", mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_postreport.html'

    script:
    """
    multiqc --filename multiqc_postreport.html .
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
    trim_galore_ch.view()
    host_remove_ch = hostremoval(trim_galore_ch)
    postfastqc_ch = postfastqc(host_remove_ch)
    postmultiqc(postfastqc_ch.collect())

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
                  
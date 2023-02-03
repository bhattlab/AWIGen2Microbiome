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


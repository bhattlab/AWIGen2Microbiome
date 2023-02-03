params.trimgalore_quality =  30
params.trimgalore_min_read_length = 60

process trimgalore {
    tag "TRIMGALORE on $sample_id"

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



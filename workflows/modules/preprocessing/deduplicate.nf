process deduplicate {
    tag "DEDUPLICATION of reads on $sample_id"

    input:
    tuple val(sample_id), path(reads)
  
    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}.fastq.gz"), emit: reads
    tuple val(sample_id), path("counts_${sample_id}.txt"), emit: stats

    script:
    """
    hts_SuperDeduper -1 ${reads[0]} -2 ${reads[1]} -f ${sample_id} -F
    count_deduplicate.py ${sample_id} stats.log > counts_${sample_id}.txt
    """
}


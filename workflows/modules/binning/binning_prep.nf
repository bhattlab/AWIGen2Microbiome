
process binning_prep {
	tag "BINNING_PREP on $sample_id"

	input:
	tuple val(sample_id), path(reads)
	tuple val(sample_id), path(contigs)

	output:
	tuple val(sample_id), path("${sample_id}.depth.txt"), emit: depth

	script:
	"""
	
	# make an index for the contigs and align the reads to get the depth
	mkdir -p ./idx
	mv $contigs ./idx
	bwa index ./idx/${contigs}
	
	bwa mem -t $task.cpus ./idx/${contigs} ${reads[0]} ${reads[1]} | samtools sort --threads $task.cpus > align.bam
	samtools sort -@ $task.cpus align.bam

	# what about singles? guess we ignore them for now?
	# seems that like jgi_summarize can take several bam files?

	jgi_summarize_bam_contig_depths --outputDepth ${sample_id}.depth.txt --pairedContigs ${sample_id}.paired.txt --minContigLength 1000 --minContigDepth 1 --percentIdentity 50 ./align.bam
	"""
}


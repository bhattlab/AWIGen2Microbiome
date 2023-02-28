process concoct {
	publishDir params.outdir + "/binning/concoct/", mode: params.publish_mode
	tag "CONCOCT on $sample_id"

	input:
	tuple val(sample_id), path(info)

	output:
	tuple val(sample_id), path("concoct_${sample_id}/"), emit: bins

	script:
	"""
	# concoct runs in multiple steps
	cut_up_fasta.py ${info[1]} -c 10000 -o 0 --merge_last \
		-b ${sample_id}_contigs_10K.bed > ${sample_id}_contigs_10K.fa
	samtools index -@ $task.cpus ${info[0]}
	concoct_coverage_table.py ${sample_id}_contigs_10K.bed \
		${info[0]} > ${sample_id}_coverage_table.tsv
	concoct --composition_file ${sample_id}_contigs_10K.fa \
		--coverage_file ${sample_id}_coverage_table.tsv \
		-b concoct_${sample_id}/ --threads $task.cpus
	merge_cutup_clustering.py concoct_${sample_id}/clustering_gt1000.csv > \
		concoct_${sample_id}/clustering_merged.csv
	extract_fasta_bins.py ${info[1]} concoct_${sample_id}/clustering_merged.csv \
		--output_path concoct_${sample_id}
	"""
}

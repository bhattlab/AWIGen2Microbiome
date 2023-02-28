process dastool {
	publishDir params.outdir + "/binning/dastool", mode: params.publish_mode
	tag "DASTool on $sample_id"

	input:
	tuple val(sample_id), path(info)

	output:
	tuple val(sample_id), path("dastool_${sample_id}"), emit: bins

	script:
	"""
	Fasta_to_Contig2Bin.sh -e fa -i ${info[0]} > metabat.tsv
	Fasta_to_Contig2Bin.sh -e fasta -i ${info[1]} > maxbin.tsv
	Fasta_to_Contig2Bin.sh -e fa -i ${info[2]} > concoct.tsv
	DAS_Tool -l metabat,maxbin,concoct --search_engine diamond \
		--threads $task.cpus --write_bins --write_unbinned \
		-i metabat.tsv,maxbin.tsv,concoct.tsv \
		-c ${info[3]} -o dastool_${sample_id} || true
	mv dastool_${sample_id}_DASTool_bins dastool_${sample_id}
	mv dastool_${sample_id}_DASTool_contig2bin.tsv ./dastool_${sample_id}
	mv dastool_${sample_id}_DASTool.log ./dastool_${sample_id}
	mv dastool_${sample_id}_DASTool_summary.tsv ./dastool_${sample_id}
	"""

}

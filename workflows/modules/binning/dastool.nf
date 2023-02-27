process dastool {
	publishDir params.outdir + "/binning/", pattern: 'dastool_*', mode: params.publish_mode
        tag "DASTool on $sample_id"

	input:
	tuple val(sample_id), path(info)

	output:
	path 'dastool_${sample_id}'

	script:
	"""
	Fasta_to_Scaffolds2Bin.sh -e fa -i ${info[0]} > metabat.tsv
        Fasta_to_Scaffolds2Bin.sh -e fasta -i ${info[1]} > maxbin.tsv
	DAS_Tool -i metabat.tsv,maxbin.tsv \
        -l ${info[0]},${info[1]} -c ${info[2]} -o dastool_${sample_id} \
        --search_engine diamond --threads $task.cpus --write_bins 1 --write_unbinned 1 || true
	"""

}

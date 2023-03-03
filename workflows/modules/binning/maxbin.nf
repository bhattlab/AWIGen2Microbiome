process maxbin {
	publishDir params.outdir + "/binning/maxbin/", mode: params.publish_mode
	tag "MAXBIN on $sample_id"

	input:
	tuple val(sample_id), path(contigs)
	tuple val(sample_id), path(depth)

	output:
	tuple val(sample_id), path("maxbin_${sample_id}/"), emit: bins

	script:
	"""
	# adjust the depth file
	mkdir -p maxbin_${sample_id}
	cut -f 1,4 ${depth} | tail -n +2 > ${depth}.adjusted
	run_MaxBin.pl -contig ${contigs} -out maxbin_${sample_id}/maxbin_bins \
		-abund ${depth}.adjusted -thread $task.cpus || true
	"""
}

process maxbin {
	tag "MAXBIN on $sample_id"

	input:
	tuple val(sample_id), path(info)

	output:
	tuple val(sample_id), path("maxbin_${sample_id}/*.fasta"), emit: bins

	script:
	"""
	# adjust the depth file
	mkdir -p maxbin_${sample_id}
	cut -f 1,4 ${info[0]} | tail -n +2 > ${info[0]}.adjusted
	run_MaxBin.pl -contig ${info[1]} -out maxbin_${sample_id}/maxbin_bins \
        	-abund ${info[0]}.adjusted -thread $task.cpus || true
        
	"""
}

params.motus_min_len_align_length = 75
params.motus_map_mgs_cutoff = 3

process motus {
	publishDir params.outdir + "/classification/motus", mode: params.publish_mode
	tag "mOTUs for $sample_id"

	input:
	tuple val(sample_id), path(reads)
	path(motus_db_path)

	output:
	path "motus_${sample_id}.out", emit: motus_res

	script:
	"""
	motus profile -n ${sample_id} -t $task.cpus -c \
		-db ${motus_db_path} \
		-l ${params.motus_min_len_align_length} -g ${params.motus_map_mgs_cutoff} \
		-f ${reads[0]} -r ${reads[1]} -s ${reads[2]} -o motus_${sample_id}.out
	"""
}

process collate_motus { 
	publishDir params.outdir + "/classification", mode: params.publish_mode

	input:
	path(motus_res)

	output:
	path "motus_all.tsv"

	script:
	"""
	mkdir -p tmp
	
	mv $motus_res ./tmp
	motus merge -o motus_all.tsv -d ./tmp
	"""
}


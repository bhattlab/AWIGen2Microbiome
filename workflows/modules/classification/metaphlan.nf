params.metaphlan_db = "/labs/asbhatt/wirbel/utils/metaphlan_db"

process metaphlan {
	publishDir params.outdir + "/classification/metaphlan", mode: params.publish_mode
	tag "MetaPhlAn for $sample_id"

	input:
	tuple val(sample_id), path(reads)

	output:
	path "metaphlan_${sample_id}.out"

	script:
	"""
	mkdir -p tmp/
	metaphlan ${reads[0]},${reads[1]},${reads[2]} \
		--input_type fastq --nrproc $task.cpus \
		--bowtie2db ${params.metaphaln_db} --tmp_dir tmp/ \
		-o metaphlan_${sample_id}.out
	"""
}



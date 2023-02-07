params.metaphlan_db_path = "/labs/asbhatt/wirbel/utils/metaphlan_db"

process metaphlan {
	publishDir params.outdir + "/classification/metaphlan", mode: params.publish_mode
	tag "MetaPhlAn for $sample_id"

	input:
	tuple val(sample_id), path(reads)
	path metaphlan_db_path

	output:
	path "metaphlan_${sample_id}.out"

	script:
	"""
	mkdir -p tmp/
	metaphlan ${reads[0]},${reads[1]},${reads[2]} \
		--input_type fastq --nproc $task.cpus \
		--bowtie2out ${sample_id}.bowtie2.bz2 \
		--bowtie2db ${metaphlan_db_path} --tmp_dir tmp/ \
		-o metaphlan_${sample_id}.out
	"""
}

process collate_metaphlan {
	publishDir params.outdir + "/classification", mode: params.publish_mode

	input:
	path(metaphlan_res)

	output:
	path "metaphlan_all.tsv"
	
	script:
	"""
	combine_metaphlan.py $metaphlan_res > metaphlan_all.tsv
	"""
}


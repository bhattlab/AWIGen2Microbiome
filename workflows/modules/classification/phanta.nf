process phanta {
	publishDir params.outdir + "/classification/phanta", mode: params.publish_mode
        tag "PHANTA for $sample_id"
	
	input:
	tuple val(sample_id), path(reads)
	path phanta_db_path

	output:
	path "phanta_${sample_id}.out", emit: phanta_res

	script:
	"""
	kraken2 --db ${phanta_db_path} --threads ${task.cpus} --output phanta_${sample_id}.krak \
		--report phanta_${sample_id}.out \
		--report-minimizer-data --paired \
		--gzip-compressed ${reads[0]} ${reads[1]} --confidence 0.1
	"""

}


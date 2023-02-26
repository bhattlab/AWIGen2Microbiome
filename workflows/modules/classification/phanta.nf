process phanta {
  publishDir params.outdir + "/classification/phanta", mode: params.publish_mode
  tag "PHANTA for $sample_id"
	
	input:
	tuple val(sample_id), path(reads)
	path phanta_db_path

	output:
	path "phanta_${sample_id}.out.*", emit: phanta_res

	script:
	"""
	kraken2 --db ${phanta_db_path} --threads ${task.cpus} \
		--output phanta_${sample_id}.krak \
		--report phanta_${sample_id}.out \
		--report-minimizer-data --paired \
		--gzip-compressed ${reads[0]} ${reads[1]} --confidence 0.1

	phanta_filter_kraken.py phanta_${sample_id}.out ${phanta_db_path} \
		${params.phanta_threshold}

	if [ ! -f phanta_${sample_id}.out.bracken ]
	then
		bracken -d {phanta_db_path} \
			-i phanta_${sample_id}.out.filtered \
    	-o phanta_${sample_id}.brack \
    	-w phanta_${sample_id}.out.bracken \
    	-r ${params.phanta_readlen} \
    	-l 'S' -t ${params.phanta_threshold}
  fi
	"""

}

process collate_phanta { 
	publishDir params.outdir + "/classification", mode: params.publish_mode

	input:
	path(phanta_res)

	output:
	path "phanta_all.tsv"

	script:
	"""
	mkdir -p bracken
	mkdir -p kraken

	mv *.out.bracken ./bracken
	mv *.out.filtered ./kraken

	combine_phanta.py ./bracken ./kraken
	"""
}


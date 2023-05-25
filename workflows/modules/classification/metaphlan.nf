params.metaphlan_db_path = "/labs/asbhatt/wirbel/utils/metaphlan_db"

process metaphlan {
	publishDir params.outdir + "/classification/metaphlan", mode: params.publish_mode, pattern: "metaphlan_*"
	tag "MetaPhlAn for $sample_id"

	input:
	tuple val(sample_id), path(reads)
	path metaphlan_db_path

	output:
	path "metaphlan_${sample_id}.out", emit: metaphlan_res
	path "versions_metaphlan_a.yml", emit: versions

	script:
	"""
	mkdir -p tmp/
	metaphlan ${reads[0]},${reads[1]},${reads[2]} \
		--input_type fastq --nproc $task.cpus \
		--bowtie2out ${sample_id}.bowtie2.bz2 \
		--bowtie2db ${metaphlan_db_path} --tmp_dir tmp/ \
		--unclassified_estimation \
		-o metaphlan_${sample_id}.out

	cat <<-END_VERSIONS > versions_metaphlan_a.yml
	"${task.process}":
	    bowtie2: \$( bowtie2 --version | head -n 1 | sed "s/.*version //g" )
	    python: \$( python --version | sed -e "s/Python //g" )
	    MetaPhlAn: \$( metaphlan --version | sed -e "s/MetaPhlAn version //g" )
	    MetaPhlAnDB: \$( cat ${metaphlan_db_path}/mpa_latest )
	END_VERSIONS
	"""
}

process collate_metaphlan {
	publishDir params.outdir + "/classification", mode: params.publish_mode, pattern: 'metaphlan_all*'

	input:
	path(metaphlan_res)

	output:
	path "metaphlan_all.tsv", emit: metaphlan_all
	path "versions_metaphlan_final.yml", emit: versions
	
	script:
	"""
	mkdir -p mpa_all
	mv $metaphlan_res ./mpa_all
	combine_metaphlan.py ./mpa_all

	cat <<-END_VERSIONS > versions_metaphlan_final.yml
	"${task.process}":
	    metaphlan_all.tsv: \$( md5sum metaphlan_all.tsv | sed -e "s/\\s.*//g" )
	END_VERSIONS

	"""
}


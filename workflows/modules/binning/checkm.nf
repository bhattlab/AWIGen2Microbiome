process checkm {
	publishDir params.outdir + "/binning/", mode: params.publish_mode
	tag "CheckM on $sample_id"

	input:
	tuple val(sample_id), path(dastool_bins)
	path(path_checkm_db)

	output:
	tuple val(sample_id), path("checkm_${sample_id}.tsv")
	path 'checkm_${sample_id}'

	shell:
	"""
	checkm data setRoot ${path_checkm_db}
	checkm lineage_wf -t $task.cpus -x fa \
		--tab_table -f checkm_${sample_id}.tsv \
		${dastool_bins} checkm_${sample_id}
	"""
}
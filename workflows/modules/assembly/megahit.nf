process megahit {
	publishDir params.outdir + "/assembly/megahit/", , mode: params.publish_mode
	tag "MEGAHIT on $sample_id"

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("megahit_${sample_id}"), emit: full_output
	tuple val(sample_id), path("megahit_${sample_id}/${sample_id}.fa"), emit: contigs

	shell:
	"""
	megahit -1 ${reads[0]} -2 ${reads[1]} -o megahit_${sample_id} --out-prefix ${sample_id}
	"""
}

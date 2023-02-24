process megahit {
	publishDir params.outdir + "/assembly/megahit_contigs/", pattern: '*.contigs.fa', mode: params.publish_mode
	tag "MEGAHIT on $sample_id"

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}.contigs.fa"), emit: contigs

	shell:
	"""
	megahit -1 ${reads[0]} -2 ${reads[1]} --out-prefix ${sample_id}
  	
	mv megahit_out/${sample_id}.contigs.fa ./
	"""
}

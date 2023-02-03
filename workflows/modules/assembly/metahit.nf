process megahit {
  publishDir params.outdir + "/assembly/${sample_id}", pattern: 'megahit_out/*contigs.fa', mode: params.publish_mode

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("megahit_out/${sample_id}.contigs.fa")

  shell:
  """
  megahit -1 ${reads[0]} -2 ${reads[1]} -r ${reads[2]} --out-prefix ${sample_id}
  """

}

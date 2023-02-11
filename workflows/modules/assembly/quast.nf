process quast {
  tag "QUAST on $sample_id"

  input:
  tuple val(sample_id), path(contigs)

  output:
  path("${sample_id}_report.tsv"), emit: quast_res

  shell:
  """
  quast.py -o quast ${contigs} --fast
  mv ./quast/report.tsv ./${sample_id}_report.tsv
  """
}

process combine_quast {
	publishDir params.outdir + "/assembly/", mode: params.publish_mode

	input:
	path(results)

	output:
	path("quast_report.tsv")
	
	shell:
	"""
	mkdir -p quast_all
	mv $results ./quast_all
	combine_quast.py ./quast_all
	"""
}


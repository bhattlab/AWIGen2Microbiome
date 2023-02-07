process quast {
  
  input:
  tuple val(sample_id), path(contigs)

  output:
  path("${sample_id}_report.tsv"), emit: quast_res

  shell:
  """
  quast.py -o quast ${contigs} --fast
  mv ./quast/transposed_report.tsv ./${sample_id}_report.tsv
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
	#!/bin/bash
	set -u
	set -e
	
	for f in $results
	do
		if [[ ! -f quast_report.tsv ]]
		then
			cat $f > quast_report.tsv
		else
			tail -n +2 $f >> quast_report.tsv
		fi
	done
	"""
}


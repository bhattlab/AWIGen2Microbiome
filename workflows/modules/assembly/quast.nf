process quast {
  input:
  tuple val(sample_id), path(contigs)

  output:
  tuple val(sample_id), path(quast)

  shell:
  """
  quast.py -o quast ${contigs} --fast
  """
}


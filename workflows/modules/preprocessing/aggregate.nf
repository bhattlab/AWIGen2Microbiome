process aggregatereports {
  publishDir params.outdir + "/stats", mode: params.publish_mode

  input:
  path '*.txt'

  output:
  path 'read_counts.tsv'

  script:
  """
  set -u
  set -e
  
  get_counts.py read_counts.tsv *.txt

  """
}


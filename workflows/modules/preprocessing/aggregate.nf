process aggregatereports {
  publishDir params.outdir + "/stats", mode: params.publish_mode

  input:
  path(stats)
  path(read_location)

  output:
  path('read_counts.tsv')
  path('count_plot.pdf')
  path('preprocessed_reads.csv')

  script:
  """
  set -u
  set -e
  
  get_counts.py read_counts.tsv $stats
  plot_counts.py read_counts.tsv

  echo "Sample ID, Forward, Reverse, Orphans" > preprocessed_reads.csv
  cat $read_location >> preprocessed_reads.csv
  """
}


process aggregatereports {
  publishDir params.outdir + "/stats", mode: params.publish_mode

  input:
  path(stats)

  output:
  path('read_counts.tsv')
  path('count_plot.pdf')

  script:
  """
  set -u
  set -e
  
  get_counts.py read_counts.tsv $stats
  plot_counts.py read_counts.tsv
  """
}


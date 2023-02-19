process aggregatereports {
  publishDir params.outdir + "/stats", mode: params.publish_mode

  input:
  path(stats)
  val(read_info)

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
  for tuple in \$(cat $read_info)
  do
    echo "\${tuple[0]},\${tuple[1]},\${tuple[2]},\${tuple[3]}" >> preprocessed_reads.csv
  done

  """
}


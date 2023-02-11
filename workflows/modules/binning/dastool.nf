process dastool {
	publishDir params.outdir + "/binning/dastool/", pattern: 'dastool_*', mode: params.publish_mode
        tag "DASTool on $sample_id"

	input:
	tuple val(sample_id) path(metabat_bins)
	tuple val(sample_id), path(maxbin_bins)
	tuple val(sample_id), path(contigs)

	output:
	path("dastool_${sample_id}")

	script:
	"""
	mkdir -p maxbin
	mv $maxbin_bins ./maxbin
	mkdir -p metabat
	mv $metabat_bins ./metabat
	Fasta_to_Scaffolds2Bin.sh -e fa -i ./metabat > metabat.tsv
        Fasta_to_Scaffolds2Bin.sh -e fasta -i ./maxbin > maxbin.tsv
	DAS_Tool -i metabat.tsv,maxbin.tsv \
        -l metabat,maxbin -c ${contigs} -o dastool_${sample_id} \
        --search_engine diamond --threads $task.cpus --write_bins 1 --write_unbinned 1 || true
	"""

}

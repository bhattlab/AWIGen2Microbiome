
process prodigal { 
	publishDir params.outdir + "/assembly", pattern: "prodigal_*", mode: params.publish_mode
	tag "PRODIGAL on $sample_id"

	input:
	tuple val(sample_id), path(contigs)

	output:
	

	shell:
	"""
	mkdir -p prodigal_${sample_id}
	prodigal -i $contigs -f gff -o prodigal_${sample_id}/${sample_id}prodigal_out.gff \
		--mrna_file prodigal_${sample_id}/${sample_id}_nucleotides.fna \
		--protein_file prodigal_${sample_id}/${sampe_id}_proteins.faa \
		-p anon
	"""
}


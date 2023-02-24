process metabat {
	tag "METABAT on $sample_id"

	input:
	tuple val(sample_id), path(info)

	output:
	tuple val(sample_id), path("metabat_${sample_id}/*fa"), emit: bins

	script:
	"""
	# actual metabat binning
	metabat2 --seed 1 -t $task.cpus \
		--unbinned \
		--inFile ${info[1]} \
		--outFile 'metabat_${sample_id}/metabat_bins' \
		--abdFile ${info[0]}
        
	# if no bins produced, copy contigs to bin.unbinned
        if [ \$(ls metabat_${sample_id} | wc -l ) == "0" ]; then
            cp ${info[1]} metabat_{sample_id}/metabat_bins.unbinned.fa
        fi

        # check for bin.tooShort.fa thats empty
        if [ -f metabat_{sample_id}/metabat_bins.tooShort.fa ]; then
            echo "Found bin.tooShort.fa"
            if [ \$(cat metabat_{sample_id}/metabat_bins.tooShort.fa | wc -l ) == "0" ]; then
                echo "Removing bin.tooShort.fa"
                rm metabat_{sample_id}/metabat_bins.tooShort.fa
            fi
        fi
	"""
}

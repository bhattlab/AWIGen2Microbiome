#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { combine_fastqs } from './modules/input/combine_fastqs'

/* PREPROCESSING
 * Processes for fastqc, multiqc, and read processing, including
 * deduplication, trimming, human read removal
*/

include { fastqc } from  './modules/preprocessing/fastqc'
include { postfastqc } from './modules/preprocessing/fastqc'
include { multiqc } from './modules/preprocessing/multiqc'
include { postmultiqc } from './modules/preprocessing/multiqc'
include { deduplicate } from './modules/preprocessing/deduplicate'
include { trimgalore } from './modules/preprocessing/trimgalore'
include { hostremoval } from './modules/preprocessing/hostremoval'
include { aggregatereports } from './modules/preprocessing/aggregate'

workflow {
	ch_sample_ids = Channel.fromList(file(params.samples).readLines())
	ch_read_pairs = combine_fastqs(ch_smaple_ids, params.read_location)
		.map { tuple(it[0], [it[1], it[2]]) }
	
	// PREPROCESSING
	ch_fastqc = fastqc(ch_read_pairs)
	multiqc(ch_fastqc.collect())
	ch_deduplicated = deduplicate(ch_read_pairs)
	ch_trim_galore = trimgalore(ch_deduplicated.reads, ch_deduplicated.stats)
	ch_host_remove = hostremoval(ch_trim_galore.reads, 
					ch_trim_galore.stats, 
					params.host_genome_location, 
					params.bwa_index_base)
	ch_postfastqc = postfastqc(ch_host_remove.reads)
	postmultiqc(ch_postfastqc.collect())
	aggregatereports(ch_host_remove.stats.collect(), ch_host_remove.read_loc.collect())
}

workflow.onComplete {
	log.info ( workflow.success ? "\nPreprocessing is done! Yay!\n" : "Oops .. something went wrong with the preprocessing" )
}

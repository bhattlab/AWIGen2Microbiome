#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
	/* FIXME
	The input should be a bit more customizable... will have to think about this in the future
	for example: what about other file endings? samples sequenced on several lanes?
	for the future, though
	*/
	Channel
		.fromFilePairs(params.reads, checkIfExists: true)
		.set { ch_read_pairs }

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

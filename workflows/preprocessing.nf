#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* PREPROCESSING
 * Processes for fastqc, multiqc, and read processing, including
 * deduplication, trimming, human read removal
*/

include { fastqc } from  "./modules/preprocessing/fastqc"
include { postfastqc } from "./modules/preprocessing/fastqc"
include { multiqc } from "./modules/preprocessing/multiqc"
include { postmultiqc } from "./modules/preprocessing/multiqc"
include { deduplicate } from "./modules/preprocessing/deduplicate"
include { trimgalore } from "./modules/preprocessing/trimgalore"
include { hostremoval } from "./modules/preprocessing/hostremoval"
include { aggregatereports } from "./modules/preprocessing/aggregate"

workflow {
	/* FIXME
	The input should be a bit more customizable... will have to think about this in the future
	for example: what about other file endings? samples sequenced on several lanes?
	for the future, though
	*/
	Channel
		.fromFilePairs(params.reads, checkIfExists: true)
		.set { read_pairs_ch }

	// PREPROCESSING
	fastqc_ch = fastqc(read_pairs_ch)
	multiqc(fastqc_ch.collect())
	deduplicated_ch = deduplicate(read_pairs_ch)
	trim_galore_ch = trimgalore(deduplicated_ch.reads, deduplicated_ch.stats)
	host_remove_ch = hostremoval(trim_galore_ch.reads, 
					trim_galore_ch.stats, 
					params.host_genome_location, 
					params.bwa_index_base)
	postfastqc_ch = postfastqc(host_remove_ch.reads)
	postmultiqc(postfastqc_ch.collect())
	aggregatereports(host_remove_ch.stats.collect())
	host_remove_ch.reads.collect().view()
}

workflow.onComplete {
	log.info ( workflow.success ? "\nPreprocessing is done! Yay!\n" : "Oops .. something went wrong with the preprocessing" )
}
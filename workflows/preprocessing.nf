#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { input_raw } from './modules/input/input_raw'

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
	ch_read_pairs = input_raw()
	ch_versions = Channel.empty()
	
	// PREPROCESSING
	// FASTQC
	ch_fastqc = fastqc(ch_read_pairs)
	ch_versions = ch_versions.mix(ch_fastqc.versions.first())

	// MULTIQC
	ch_multiqc_pre = multiqc(ch_fastqc.prefastqc.collect())
	ch_versions = ch_versions.mix(ch_multiqc_pre.versions.first())

	// DEDUPLICATION
	ch_deduplicated = deduplicate(ch_read_pairs)
	ch_versions = ch_versions.mix(ch_deduplicated.versions.first())

	// TRIMMING
	ch_trim_galore = trimgalore(ch_deduplicated.reads, ch_deduplicated.stats)
	ch_versions = ch_versions.mix(ch_trim_galore.versions.first())

	// HOST REMOVAL
	ch_host_remove = hostremoval(ch_trim_galore.reads, 
					ch_trim_galore.stats, 
					params.host_genome_location, 
					params.bwa_index_base)
	ch_versions = ch_versions.mix(ch_host_remove.versions.first())

	// FASTQC again
	ch_postfastqc = postfastqc(ch_host_remove.reads)

	// MULTIQC again
	postmultiqc(ch_postfastqc.postfastqc.collect())

	// AGGREGATE REPORTS
	aggregatereports(ch_host_remove.stats.collect(), 
					ch_host_remove.read_loc.collect())

	// VERSION output
	ch_versions
		.unique()
		.collectFile(name: params.outdir + 'versions_preprocessing.yml')

}

workflow.onComplete {
	log.info ( workflow.success ? "\nPreprocessing is done! Yay!\n" : "Oops .. something went wrong with the preprocessing" )
}

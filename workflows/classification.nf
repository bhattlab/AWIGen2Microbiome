#!/usr/bin/env nextflow

nextflow.enable.dsl=2


include { input_check } from './modules/input/input_check'

/* CLASSIFICATION
 * Processes for motus, metaphlan, phanta (including kraken2)
*/

include { motus } from  './modules/classification/motus'
include { collate_motus } from './modules/classification/motus'
include { metaphlan } from './modules/classification/metaphlan'
include { collate_metaphlan } from './modules/classification/metaphlan'
include { phanta } from './modules/classification/phanta'

workflow {

	// Input
	ch_processed_reads = input_check()
	// CLASSIFICATION
	if (params.run_motus) {
		ch_motus = motus(ch_processed_reads.reads, params.motus_db_path)
		ch_all_motus = collate_motus(ch_motus.motus_res.collect(), params.motus_db_path)
	}
	if (params.run_metaphlan) {
		ch_metaphlan = metaphlan(ch_processed_reads.reads, params.metaphlan_db_path)
		ch_all_metaphlan = collate_metaphlan(ch_metaphlan.metaphlan_res.collect())
	}
	if (params.run_phanta) {
		ch_phanta = phanta(ch_processed_reads.reads, params.phanta_db_path)
	//	ch_all_phanta = collate_phanta(ch_phanta.phanta_res.collect())
	}
}

workflow.onComplete {
	log.info ( workflow.success ? "\nClassification is done! Yay!\n" : "Oops .. something went wrong" )
}

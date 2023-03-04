#!/usr/bin/env nextflow

nextflow.enable.dsl=2


include { input_check } from './modules/input/input_check'
include { input_assembly } from './modules/input/input_assembly'

/* ASSEMBLY & BINNING
 * use a complete workflow so that you don't have to wait for 
 * all samples to finish each step before channels can be combined
*/

include { assembly_binning } from './subworkflows/assembly_binning'

workflow {

	// Input
	ch_processed_reads = input_check()
	// ASSEMBLY
	ch_bin = assembly_binning(ch_processed_reads.reads)

	// DAStool
	ch_dastool = dastool(
		ch_bin.contigs,
		ch_bin.metabat, 
		ch_bin.maxbin, 
		ch_bins.concoct)

	ch_checkm = checkm(ch_dastool.bins, params.checkm_db_path)
	// ch_quast_res = combine_quast(ch_quast.quast_res.collect())
}

workflow.onComplete {
	log.info ( workflow.success ? "\nAssembly and binning is done! Yay!\n" : "Oops .. something went wrong" )
}


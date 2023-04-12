#!/usr/bin/env nextflow

nextflow.enable.dsl=2


include { input_check } from './modules/input/input_check'

/* ASSEMBLY
 run megahit, prodigal, and quast
*/

include { megahit } from './modules/assembly/megahit'

workflow {

	// Input
	ch_processed_reads = input_check()
	// ASSEMBLY
	ch_megahit = megahit(ch_processed_reads)
	ch_quast = quast(ch_megahit.contigs)
	ch_prodigal = prodigal(ch_megahit.contigs)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nAssembly and binning is done! Yay!\n" : "Oops .. something went wrong" )
}

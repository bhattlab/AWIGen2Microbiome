#!/usr/bin/env nextflow

nextflow.enable.dsl=2


include { input_check } from './modules/input/input_check'

/* ASSEMBLY
 * Runs megahit on paired and orphan reads, then uses QUAST to measure
 * assembly quality. Combines quast reports into a single report.
*/

include { megahit } from './modules/assembly/megahit'
include { prodigal } from './modules/assembly/prodigal'
include { quast } from './modules/assembly/quast'
include { combine_quast } from './modules/assembly/quast'

/* BINNING
 * Runs metabat2 and maxbin on the contigs, then DAS tool (for later)
 * Check quality by CheckM or sth (maybe others)?
*/

include { binning_prep } from './modules/binning/binning_prep'
include { metabat } from './modules/binning/metabat'
include { maxbin } from './modules/binning/maxbin'
include { dastool } from './modules/binning/dastool'


workflow {

	// Input
	ch_processed_reads = input_check()
	// ASSEMBLY
	megahit_ch = megahit(host_remove_ch.reads)
	quast_ch = quast(megahit_ch.contigs)
	prodigal_ch = prodigal(megahit_ch.contigs)
	quast_res_ch = combine_quast(quast_ch.quast_res.collect())
    
	// BINNING
	binning_prep_ch = binning_prep(host_remove_ch.reads, megahit_ch.contigs)
	metabat_bins_ch = metabat(megahit_ch.contigs, binning_prep_ch.depth)
	maxbin_bins_ch = maxbin(megahit_ch.contigs, binning_prep_ch.depth)
	dastool_ch = dastool(metabat_bins_ch.bins, maxbin_bins_ch.bins, megahit_ch.contigs)

	// checkm
	// other things?
}

workflow.onComplete {
	log.info ( workflow.success ? "\nAssembly and binning is done! Yay!\n" : "Oops .. something went wrong" )
}


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
	ch_megahit = megahit(ch_processed_reads.reads)
	ch_quast = quast(ch_megahit.contigs)
	ch_prodigal = prodigal(ch_megahit.contigs)
	ch_quast_res = combine_quast(ch_quast.quast_res.collect())
    
	// BINNING
	ch_binning_prep = binning_prep(ch_processed_reads.reads, ch_megahit.contigs)
	ch_metabat_bins = metabat(ch_megahit.contigs, ch_binning_prep.depth)
	ch_maxbin_bins = maxbin(ch_megahit.contigs, ch_binning_prep.depth)
	ch_dastool = dastool(ch_metabat_bins.bins, ch_maxbin_bins.bins, ch_megahit.contigs)

	// checkm
	// other things?
}

workflow.onComplete {
	log.info ( workflow.success ? "\nAssembly and binning is done! Yay!\n" : "Oops .. something went wrong" )
}


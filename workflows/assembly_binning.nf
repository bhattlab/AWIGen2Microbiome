#!/usr/bin/env nextflow

nextflow.enable.dsl=2


include { input_check } from './modules/input/input_check'
include { input_assembly } from './modules/input/input_assembly'

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
include { concoct } from './modules/binning/concoct'
include { dastool } from './modules/binning/dastool'
include { checkm } from './modules/binning/checkm'

workflow {

	// Input
	ch_processed_reads = input_check()
	// ASSEMBLY
	// ch_megahit = megahit(ch_processed_reads.reads)
	ch_megahit = input_assembly()
	ch_quast = quast(ch_megahit.contigs)
	ch_prodigal = prodigal(ch_megahit.contigs)
	ch_quast_res = combine_quast(ch_quast.quast_res.collect())
    
	// BINNING
	ch_binning_prep = binning_prep(ch_processed_reads.reads, ch_megahit.contigs)
	ch_binning = ch_binning_prep.depth
		.concat(ch_megahit.contigs)
		.groupTuple()
	
	ch_metabat_bins = metabat(ch_binning)
	ch_maxbin_bins = maxbin(ch_binning)
	ch_binning_concoct = ch_binning_prep.idx
		.concat(ch_megahit.contigs)
		.groupTuple()
	ch_concoct_bins = concoct(ch_binning_concoct)
	
	// combine the metabat, maxbin, contigs channels
	ch_final = ch_metabat_bins.bins
		.concat(ch_maxbin_bins.bins)
		.concat(ch_concoct_bins.bins)
		.concat(ch_megahit.contigs).groupTuple()
	ch_dastool = dastool(ch_final)

	// checkm
	ch_checkm = checkm(ch_dastool.bins, params.checkm_db_path)
	// other things?
}

workflow.onComplete {
	log.info ( workflow.success ? "\nAssembly and binning is done! Yay!\n" : "Oops .. something went wrong" )
}


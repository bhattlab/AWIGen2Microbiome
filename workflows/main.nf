#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// params.reads = "/labs/asbhatt/dmaghini/projects/awigen2/workflow_planning/00.test_samples/*_{1,2}.fq.gz"
// params.reads = "/labs/asbhatt/dmaghini/tools/nextflow_tutorial/indata/*{1,2}.fq.gz"
params.reads = "./test_dataset/raw_data/*{1,2}.fq.gz"

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

/* CLASSIFICATION
 * Processes for motus, metaphlan, phanta (including kraken2)
*/

include { motus } from  "./modules/classification/motus"
include { collate_motus } from "./modules/classification/motus"
include { metaphlan } from "./modules/classification/metaphlan"
include { collate_metaphlan } from "./modules/classification/metaphlan"
include { phanta } from "./modules/classification/phanta"

/* ASSEMBLY
 * Runs megahit on paired and orphan reads, then uses QUAST to measure
 * assembly quality. Combines quast reports into a single report.
*/

include { megahit } from "./modules/assembly/megahit"
include { prodigal } from "./modules/assembly/prodigal"
include { quast } from "./modules/assembly/quast"
include { combine_quast } from "./modules/assembly/quast"

/* BINNING
 * Runs metabat2 and maxbin on the contigs, then DAS tool (for later)
 * Check quality by CheckM or sth (maybe others)?
*/

include { binning_prep } from "./modules/binning/binning_prep"
include { metabat } from "./modules/binning/metabat"
include { maxbin } from "./modules/binning/maxbin"
include { dastool } from "./modules/binning/dastool"

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
    
	// CLASSIFICATION
	if (params.run_motus) {
		motus_ch = motus(host_remove_ch.reads)
		motus_all = collate_motus(motus_ch.motus_res.collect())
	}
	if (params.run_metaphlan) {
		metaphlan_ch = metaphlan(host_remove_ch.reads, params.metaphlan_db_path)
		metaphlan_all = collate_metaphlan(metaphlan_ch.metaphlan_res.collect())
	}
	if (params.run_phanta) {
		phanta_ch = phanta(host_remove_ch.reads)
	//	phanta_all = collate_phanta(phanta_ch.phanta_res.collect())
	}

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
	log.info ( workflow.success ? "\nDone! Yay!\n" : "Oops .. something went wrong" )
}

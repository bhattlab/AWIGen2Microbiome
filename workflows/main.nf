#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
 * pipeline input parameters -- This should go into a config file before long
 */
params.reads = "/labs/asbhatt/dmaghini/tools/nextflow_tutorial/indata/*{1,2}.fq.gz"


log.info """\
    METAGENOMIC PREPROCESSING - NF PIPELINE
    ===================================
    raw.reads        : ${params.reads}
    host.genome      : ${params.bwa_index_base}
    
    outdir           : ${params.outdir}
    """
    .stripIndent()

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


process aggregatereports {
  publishDir params.outdir + "/stats/read_counts", mode: params.publish_mode
  
  input:
  path rawstats
  path dedupstats
  path trimstats
  path hoststats

  output:
  path 'readcounts.tsv'

  script:
  """
  cat ${rawstats}  ${dedupstats} ${trimstats} ${hoststats} > readcounts.tsv
  """
}

/* CLASSIFICATION
 * Processes for motus, metaphlan, phanta (including kraken2)
*/

include { motus } from  "./modules/classification/motus"

/* ASSEMBLY
 * Runs megahit on paired and orphan reads, then uses QUAST to measure
 * assembly quality. Combines quast reports into a single report.
*/

include { megahit } from "./modules/assembly/megahit"
include { quast } from "./modules/assembly/quast"

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
    multiqc(fastqc_ch.logs.collect())
    deduplicated_ch = deduplicate(read_pairs_ch)
    trim_galore_ch = trimgalore(deduplicated_ch.dedupreads)
    host_remove_ch = hostremoval(trim_galore_ch.trimreads, params.host_genome_location, params.bwa_index_base)
    postfastqc_ch = postfastqc(host_remove_ch.hostremreads)
    postmultiqc(postfastqc_ch.collect())
    //aggregatereports(fastqc_ch.stats.collect(), 
    //                 deduplicated_ch.dedupstats.collect(), 
    //                 trim_galore_ch.trimstats.collect(), 
    //                 host_remove_ch.hoststats.collect())
    
    // CLASSIFICATION
    motus_ch = motus(host_remove_ch.hostremreads)
    // ASSEMBLY
    megahit_ch = megahit(host_remove_ch.hostremreads)
    quast(megahit_ch)
    
    // BINNING

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Yay!\n" : "Oops .. something went wrong" )
}

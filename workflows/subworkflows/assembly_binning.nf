
/* ASSEMBLY
 * Runs megahit on paired and orphan reads, then uses QUAST to measure
 * assembly quality. Combines quast reports into a single report.
*/

include { megahit } from '../modules/assembly/megahit'
include { prodigal } from '../modules/assembly/prodigal'
include { quast } from '../modules/assembly/quast'
include { combine_quast } from '../modules/assembly/quast'

/* BINNING
 * Runs metabat2 and maxbin on the contigs, then DAS tool (for later)
 * Check quality by CheckM or sth (maybe others)?
*/

include { binning_prep } from '../modules/binning/binning_prep'
include { metabat } from '../modules/binning/metabat'
include { maxbin } from '../modules/binning/maxbin'
include { concoct } from '../modules/binning/concoct'
include { dastool } from '../modules/binning/dastool'
include { checkm } from '../modules/binning/checkm'

workflow assembly_binning {

    take:
        input //channel: tuple: val(sample_id), path(reads)

    main:

    ch_megahit = megahit(input)
    ch_quast = quast(ch_megahit.contigs)
    ch_prodigal = prodigal(ch_megahit.contigs)

    ch_binning_prep = binning_prep(input, ch_megahit.contigs)
    ch_binning = ch_binning_prep.depth
        .concat(ch_megahit.contigs)
        .groupTuple()
    
    ch_metabat_bins = metabat(ch_binning)
    ch_maxbin_bins = maxbin(ch_binning)
    ch_binning_concoct = ch_binning_prep.bam
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
}
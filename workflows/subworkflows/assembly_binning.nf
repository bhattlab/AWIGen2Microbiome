
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

    ch_binning_prep = binning_prep(ch_megahit.reads, ch_megahit.contigs)

    ch_metabat_bins = metabat(ch_binning_prep.contigs, ch_binning_prep.depth)
    ch_maxbin_bins = maxbin(ch_binning_prep.contigs, ch_binning_prep.depth)
    ch_concoct_bins = concoct(ch_binning_prep.contigs, ch_binning_prep.bam)

    emit:
    contigs = ch_megahit.contigs
    metabat = ch_metabat_bins.bins
    maxbin = ch_maxbin_bins.bins
    concoct = ch_concoct_bins.bins


}

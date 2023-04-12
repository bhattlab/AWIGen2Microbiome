process multiqc {
    publishDir params.outdir + "/stats/multiqc_pre", mode: params.publish_mode
    tag "MULTIQC before anything"

    input:
    path '*'

    output:
    path 'multiqc_pre_report*'

    script:
    """
    multiqc --filename multiqc_pre_report.html .
    """
}

process postmultiqc {
    publishDir params.outdir + "/stats/multiqc_post", mode: params.publish_mode

    input:
    path '*'

    output:
    path 'multiqc_postreport*'

    script:
    """
    multiqc --filename multiqc_postreport.html .
    """
}


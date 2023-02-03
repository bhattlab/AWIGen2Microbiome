process multiqc {
    publishDir params.outdir + "/stats", mode: params.publish_mode
    tag "MULTIQC before anything"

    input:
    path '*'

    output:
    path 'multiqc_pre_report.html'

    script:
    """
    multiqc --filename multiqc_pre_report.html .
    """
}

process postmultiqc {
    publishDir params.outdir + "/stats", mode: params.publish_mode

    input:
    path '*'

    output:
    path 'multiqc_postreport.html'

    script:
    """
    multiqc --filename multiqc_postreport.html .
    """
}


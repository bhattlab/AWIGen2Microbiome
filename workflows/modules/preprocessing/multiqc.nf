process multiqc {
    publishDir params.outdir + "/stats/pre_multiqc/", mode: params.publish_mode, pattern: "multiqc*"
    tag "MULTIQC before anything"

    input:
    path '*'

    output:
    path 'multiqc_pre_report*', emit: premultiqc
    path "versions.yml", emit: versions

    script:
    """
    multiqc --filename multiqc_pre_report.html .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}

process postmultiqc {
    publishDir params.outdir + "/stats/post_multiqc", mode: params.publish_mode

    input:
    path '*'

    output:
    path 'multiqc_postreport*', emit: postmultiqc

    script:
    """
    multiqc --filename multiqc_postreport.html .
    """
}


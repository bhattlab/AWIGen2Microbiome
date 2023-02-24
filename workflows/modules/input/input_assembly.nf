// workflow to read a CSV file with sample names and read files
// and emit them as tuples for downstream analysis

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

workflow input_assembly {
    main:
    if(hasExtension(params.input_assembly, "csv")){
        // extracts read files from samplesheet CSV and distribute into channels
        ch_input = Channel
            .from(file(params.input_assembly))
            .splitCsv(header: true)
            .map { row ->
                    if (row.size() == 2) {
                        def sampleid = row.sampleID
                        def contigs = row.contigs ? file(row.contigs, checkIfExists: true) : false
                        // Check if given combination is valid
                        if (!contigs) exit 1, "Invalid input samplesheet: Contigs can not be empty."
                        return [ sampleid, contigs ]
                    } else {
                        exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects 2."
                    }
             }
        ch_contigs = ch_input
            .map { sampleid, contigs ->
                        return [ sampleid, contigs ]
                }
    } else {
        exit 1, "Input samplesheet should be a csv file organised like this:\n\nsampleID,contigs"
    }
   
    emit:
    contigs = ch_contigs
}

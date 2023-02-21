// workflow to read a CSV file with sample names and read files
// and emit them as tuples for downstream analysis

workflow input_check {
    main:
    if(hasExtension(params.input, "csv")){
        // extracts read files from samplesheet CSV and distribute into channels
        ch_input = Channel
            .from(file(params.input))
            .splitCsv(header: true)
            .map { row ->
                    if (row.size() == 4) {
                        def sampleid = row.sample
                        def forward = row.forward ? file(row.forward, checkIfExists: true) : false
                        def reverse = row.reverse ? file(row.reverse, checkIfExists: true) : false
                        def orphans = row.orphans ? file(row.orphans, checkIfExists: true) : false
                        // Check if given combination is valid
                        if (!forward) exit 1, "Invalid input samplesheet: Forward_reads can not be empty."
                        if (!reverse) exit 1, "Invalid input samplesheet: Reverse_reads can not be empty."
                        if (!orphans) exit 1, "Invalid input samplesheet: Orphan_reads can not be empty."
                        return [ sampleid, forward, reverse, orphans ]
                    } else {
                        exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects 4."
                    }
             }
        ch_reads = ch_input
            .map { sampleid, forward, reverse, orphans ->
                        def meta = [:]
                        return [ sampleid, [ sr1, sr2 ] ]
                }
        emit:
        reads = ch_reads
    } else {
        exit 1, "Input samplesheet should be a csv file organised like this:\n\nSample_ID,Forward_reads,Reverse_reads,Orphan_reads"
    }
}
include { RMARKDOWNNOTEBOOK as CREATE_REPORT    } from '../../modules/nf-core/rmarkdownnotebook/main'

workflow POSTPROCESS {

    take:
    ch_vcf
    ch_fst
    ch_ref
    ch_sync
    ch_split_sync
    ch_frequency
    ch_fisher

    main:

    ch_versions = Channel.empty()

    // TODO: figure out what happens if there's no fisher test

    // collapse pairwise fisher tests into single file
    ch_fisher_collapsed = ch_fisher
        .map { it[1] }
        .collectFile( name: 'fisher_all.tsv', keepHeader: true, storeDir: 'output/fishertest' )

    // build rmarkdown report
    def report_file_names = [ 'fst_file', 'fisher' ]
    ch_report = ch_vcf.map { [ it[0], file("${projectDir}/assets/assesspool_report.Rmd") ] }
    ch_input_files = ch_fst
        .combine( ch_fisher_collapsed )
        .collect()
    ch_params = ch_input_files.map {
        [report_file_names, it.collect{ it.name }].transpose().collectEntries() +
            [ 'viz_filter': params.visualize_filters ] +
            [ 'tz': TimeZone.getDefault().getID() ] // this is required because something about the container goes nuts without it
    }

    CREATE_REPORT( ch_report, ch_params, ch_input_files )

    emit:
    // // TODO nf-core: edit emitted channels

    versions = ch_versions                     // channel: [ versions.yml ]
}


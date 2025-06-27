include { RMARKDOWNNOTEBOOK as CREATE_REPORT    } from '../../modules/nf-core/rmarkdownnotebook/main'

workflow POSTPROCESS {

    take:
    ch_vcf
    ch_fst
    ch_ref
    ch_sync
    ch_split_sync
    ch_frequency

    main:

    ch_versions = Channel.empty()

    // build rmarkdown report
    def report_file_names = [ 'fst_file' ]
    ch_report = ch_vcf.map { [ it[0], file("${projectDir}/assets/assesspool_report.Rmd") ] }
    ch_input_files = ch_fst.collect()
    ch_params = ch_input_files.map {
        [report_file_names, it.collect{ it.name }].transpose().collectEntries() +
            [ 'debug': true] +
            [ 'viz_filter': params.visualizations ] +
            [ 'tz': TimeZone.getDefault().getID() ] // this is required because something about the container goes nuts without it
    }

    CREATE_REPORT( ch_report, ch_params, ch_input_files )



    emit:
    // // TODO nf-core: edit emitted channels

    versions = ch_versions                     // channel: [ versions.yml ]
}


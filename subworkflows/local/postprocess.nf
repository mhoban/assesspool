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
    ch_filter

    main:

    ch_versions = Channel.empty()

    // collapse pairwise fisher tests into single file
    ch_fisher_collapsed = ch_fisher
        .map { it[1] }
        .collectFile( name: 'fisher_all.tsv', keepHeader: true, storeDir: 'output/fishertest' )

    // build rmarkdown report
    // def file_keys = [ 'fst_file', 'fisher', 'filter' ]
    ch_report = ch_vcf.map { [ it[0], file("${projectDir}/assets/assesspool_report.Rmd") ] }
    ch_input_files = ch_fst.ifEmpty{ [] }
        .combine( ch_fisher_collapsed.ifEmpty{ [] } )
        .combine( ch_filter.ifEmpty{ [] } )
        .collect()

    ch_fst = ch_fst
        .map{ [ fst_file: it.name ] }
        .ifEmpty{[:]}


    ch_fisher_collapsed = ch_fisher_collapsed
        .map{ [ fisher: it.name ] }
        .ifEmpty{[:]}

    ch_filter = ch_filter
        .map{ [ filter: it.name ] }
        .ifEmpty{[:]}

    ch_params = ch_fst
        .mix(ch_fisher_collapsed)
        .mix(ch_filter)
        .reduce{a,b -> a+b }
        .map{ it + [ 'viz_filter': params.visualize_filters, 'tz': TimeZone.getDefault().getID() ] }

    CREATE_REPORT( ch_report, ch_params, ch_input_files )

    emit:
    // // TODO nf-core: edit emitted channels

    versions = ch_versions                     // channel: [ versions.yml ]
}


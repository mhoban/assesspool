include { RMARKDOWNNOTEBOOK as CREATE_REPORT } from '../../modules/nf-core/rmarkdownnotebook/main'
include { RIPGREP as COUNT_SNPS_FINAL        } from '../../modules/nf-core/ripgrep/main'
include { EXTRACT_SEQUENCES                   } from '../../modules/local/extractsequences/main'

workflow POSTPROCESS {

    take:
    ch_input
    ch_unfiltered
    ch_fst
    ch_ref
    ch_sync
    ch_split_sync
    ch_frequency
    ch_fisher
    ch_filter

    main:

    ch_versions = Channel.empty()

    // extract reference contigs with "strongly differentiated" snps
    if (params.extract_sequences) {
        EXTRACT_SEQUENCES( ch_ref, ch_fst, params.fst_cutoff )
        ch_versions = ch_versions.mix(EXTRACT_SEQUENCES.out.versions.first())
    }


    // collapse pairwise fisher tests into single files
    ch_fisher_collapsed = ch_fisher
        .collectFile( keepHeader: true, sort: true, storeDir: 'output/fishertest'  ){ meta, fish -> [ "${meta.id}.fisher", fish ] }
        .map{ [ it.baseName, it ] }

    // join them back to meta tags
    ch_fisher_collapsed = ch_input
        .map{ it[0] }
        .unique()
        .map{ [ it.id, it ] }
        .join( ch_fisher_collapsed )
        .map{ id, meta, f -> [ meta, f ] }

    ch_count = ch_input
        .map{ meta, input, index -> [ meta + [filter: 'cumulative'], input ] }
        .mix( ch_unfiltered.map{ meta, input, index -> [ meta + [filter: 'before'], input ] } )

    // count final filtered SNPs into map
    COUNT_SNPS_FINAL( ch_count, '^#', false )
    ch_filter_final = COUNT_SNPS_FINAL.out.txt.map{ meta, count -> [ meta, [ filter: meta.filter, count: count.text.trim() ] ] }

    // collect final filter summary into tsv files
    ch_filter_final = ch_filter_final
        .collectFile(newLine: true, sort: true ) { meta, filter -> [ "${meta.id}.final_filter", "${filter.filter}\t${filter.count}" ] }
        .map{ [ it.baseName, it ] }

    // join them back to the meta tags
    ch_filter_final = ch_input
        .map{ it[0] }
        .unique()
        .map{ [ it.id, it ] }
        .join( ch_filter_final )
        .map{ id, meta, f -> [ meta, f ] }

    // build rmarkdown report input and params

    // get report file channel as [ meta, reportfile ]
    ch_report = ch_input.map { [ it[0].id, it[0], file("${projectDir}/assets/assesspool_report.Rmd") ] }

    // generate input files channel
    ch_input_files = ch_fst.ifEmpty{ [] }
        .mix( ch_fisher_collapsed.ifEmpty{ [] } )
        .mix( ch_filter.ifEmpty{ [] } )
        .mix( ch_filter_final.ifEmpty{ [] } )
        .groupTuple()
        .map{ [it[0].id] + it[1..-1] }


    // subset the params object because there's at least one value that changes
    // every time, which invalidates the caching
    nf_params = [
        'coverage_cutoff_step',
        'max_coverage_cutoff',
        'min_coverage_cutoff',
        'visualize_filters'
    ]

    ch_params = ch_input_files.
        map{ meta, files -> [ meta, files.collect{ [ it.extension, it.name ] }.collectEntries() ] }
        .map{ meta, p -> [ meta, p + [
            nf: params.subMap(nf_params),
            tz: TimeZone.getDefault().getID()
        ]]}

    // join everything together
    ch_report
        .join( ch_params )
        .join( ch_input_files )
        .map { it[1..-1] }
        .set{ ch_report }

    // tuple val(meta), path(notebook)
    // val parameters
    // path input_files
    CREATE_REPORT(
        ch_report.map{ meta, report, params, files -> [meta, report] },
        ch_report.map{ meta, report, params, files -> params },
        ch_report.map{ meta, report, params, files -> files }
    )
    ch_versions = ch_versions.mix(CREATE_REPORT.out.versions.first())

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}


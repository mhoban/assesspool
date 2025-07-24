include { RMARKDOWNNOTEBOOK as CREATE_REPORT    } from '../../modules/nf-core/rmarkdownnotebook/main'
include { RIPGREP as COUNT_SNPS_FINAL           } from '../../modules/nf-core/ripgrep/main'

workflow POSTPROCESS {

    take:
    ch_vcf
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

    // collapse pairwise fisher tests into single files
    ch_fisher_collapsed = ch_fisher
        .collectFile( keepHeader: true, sort: false, storeDir: 'output/fishertest' ){ meta, fish -> [ "${meta.id}.fisher", fish ] }
        .map{ [ it.baseName, it ] }
    // join them back to meta tags
    ch_fisher_collapsed = ch_vcf
        .map{ it[0] }
        .unique()
        .map{ [ it.id, it ] }
        .join( ch_fisher_collapsed )
        .map{ id, meta, f -> [ meta, f ] }

    ch_count = ch_vcf
        .map{ meta, vcf, index -> [ meta + [filter: 'cumulative'], vcf ] }
        .mix( ch_unfiltered.map{ meta, vcf, index -> [ meta + [filter: 'before'], vcf ] } )

    // count final filtered SNPs into map
    COUNT_SNPS_FINAL( ch_count, '^#', false )
    ch_filter_final = COUNT_SNPS_FINAL.out.txt.map{ meta, count -> [ meta, [ filter: meta.filter, count: count.text.trim() ] ] }

    // collect final filter summary into tsv files
    ch_filter_final = ch_filter_final
        .map{ meta, count -> meta.subMap('id') }
        .unique()
        .map{ meta -> [ meta, [ filter: 'filter', count: 'count' ] ] }
        .concat( ch_filter_final )
        .collectFile(newLine: true, sort: false) { meta, filter -> [ "${meta.id}.final_filter", "${filter.filter}\t${filter.count}" ] }
        .map{ [ it.baseName, it ] }
    // join them back to the meta tags
    ch_filter_final = ch_vcf
        .map{ it[0] }
        .unique()
        .map{ [ it.id, it ] }
        .join( ch_filter_final )
        .map{ id, meta, f -> [ meta, f ] }

    // build rmarkdown report input and params

    // get report file channel as [ meta, reportfile ]
    ch_report = ch_vcf.map { [ it[0], file("${projectDir}/assets/assesspool_report.Rmd") ] }

    // generate input files channel
    ch_input_files = ch_fst.ifEmpty{ [] }
        .mix( ch_fisher_collapsed.ifEmpty{ [] } )
        .mix( ch_filter.ifEmpty{ [] } )
        .mix( ch_filter_final.ifEmpty{ [] } )
        .groupTuple()

    ch_params = ch_input_files.
        map{ meta, files -> [ meta, files.collect{ [ it.extension, it.name ] }.collectEntries() ] }
        .map{ meta, p -> [ meta, p + [
            nf: params,
            tz: TimeZone.getDefault().getID()
        ]]}

    CREATE_REPORT( ch_report, ch_params.map{meta, p -> p}, ch_input_files.map{meta, f -> f} )
    ch_versions = ch_versions.mix(CREATE_REPORT.out.versions.first())

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}


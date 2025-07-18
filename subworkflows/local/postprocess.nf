include { RMARKDOWNNOTEBOOK as CREATE_REPORT    } from '../../modules/nf-core/rmarkdownnotebook/main'
include { RIPGREP as COUNT_SNPS_FINAL           } from '../../modules/nf-core/ripgrep/main'

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

    // count final filtered SNPs into map
    COUNT_SNPS_FINAL( ch_vcf.map{ meta, vcf, index -> [ meta, vcf ] }, '^#', false )
    ch_filter_final = COUNT_SNPS_FINAL.out.txt.map{ meta, count -> [ meta, [ filter: 'cumulative', count: count.text.trim() ] ] }

    // collect final filter summary into tsv file
    ch_filter_final = ch_filter_final
        .map{ meta, count -> meta.subMap('id') }
        .unique()
        .map{ meta -> [ meta, [ filter: 'filter', count: 'count' ] ] }
        .concat( ch_filter_final )
        .collectFile(newLine: true, sort: false) { meta, filter -> [ "${meta.id}_final_filter.tsv", "${filter.filter}\t${filter.count}" ] }

    // build rmarkdown report
    // def file_keys = [ 'fst_file', 'fisher', 'filter' ]
    ch_report = ch_vcf.map { [ it[0], file("${projectDir}/assets/assesspool_report.Rmd") ] }
    ch_input_files = ch_fst.ifEmpty{ [] }
        .combine( ch_fisher_collapsed.ifEmpty{ [] } )
        .combine( ch_filter.ifEmpty{ [] } )
        .combine( ch_filter_final.ifEmpty{ [] } )
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

    ch_filter_final = ch_filter_final
        .map{ [ final_filter: it.name ] }
        .ifEmpty{[:]}

    ch_params = ch_fst
        .mix(ch_fisher_collapsed)
        .mix(ch_filter)
        .mix(ch_filter_final)
        .reduce{a,b -> a+b }
        .map{ it + [
            'nextflow_params': params,
            'viz_filter': params.visualize_filters,
            'tz': TimeZone.getDefault().getID()
        ] }

    CREATE_REPORT( ch_report, ch_params, ch_input_files )

    emit:
    // // TODO nf-core: edit emitted channels

    versions = ch_versions                     // channel: [ versions.yml ]
}


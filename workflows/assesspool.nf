/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_assesspool_pipeline'
include { FILTER                 } from '../subworkflows/local/filter.nf'
include { FILTER_SIM             } from '../subworkflows/local/filter_sim/main.nf'
include { PREPROCESS             } from '../subworkflows/local/preprocess.nf'
include { POOLSTATS              } from '../subworkflows/local/poolstats.nf'
include { POSTPROCESS             } from '../subworkflows/local/postprocess.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSESSPOOL {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    // run VCF preprocessing
    PREPROCESS( ch_samplesheet )
    ch_vcf = PREPROCESS.out.vcf
    ch_ref = PREPROCESS.out.ref
    ch_versions = ch_versions.mix(PREPROCESS.out.versions.first())

    // prepare channels
    ch_fst = Channel.empty()
    ch_filter_sim = Channel.empty()
    ch_sync = Channel.empty()
    ch_split_sync = Channel.empty()
    ch_freq = Channel.empty()
    ch_fisher = Channel.empty()
    ch_filtered = Channel.empty()

    // perform stepwise filtration
    // for visualization and evaluation


    if (params.visualize_filters) {
        FILTER_SIM( ch_vcf )
        ch_filter_sim = FILTER_SIM.out.filter_summary
    }

    // run downstream processing unless we're only
    // visualizing filters
    if (!params.filter_only) {
        // run filtering if desired
        FILTER( ch_vcf )
        ch_filtered = FILTER.out.vcf
        ch_versions = ch_versions.mix(FILTER.out.versions.first())

        // calculate pool statistics (fst, fisher, etc.)
        POOLSTATS(ch_filtered,ch_ref)
        ch_versions = ch_versions.mix(POOLSTATS.out.versions.first())

        ch_fst = POOLSTATS.out.fst
        ch_sync = POOLSTATS.out.sync
        ch_split_sync = POOLSTATS.out.split_sync
        ch_freq = POOLSTATS.out.frequency
        ch_fisher = POOLSTATS.out.fisher
    }

    // run post-processing steps
    POSTPROCESS(
        ch_filtered,
        ch_vcf,
        ch_fst,
        ch_ref,
        ch_sync,
        ch_split_sync,
        ch_freq,
        ch_fisher,
        ch_filter_sim
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'assesspool_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

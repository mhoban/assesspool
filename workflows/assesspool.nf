/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_assesspool_pipeline'
include { FILTER                 } from '../subworkflows/local/filter.nf'
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
    PREPROCESS(ch_samplesheet)
    ch_vcf = PREPROCESS.out.vcf
    ch_ref = PREPROCESS.out.ref
    ch_versions = ch_versions.mix(PREPROCESS.out.versions.first())

    // run filtering if desired
    FILTER(ch_vcf)
    ch_filtered = FILTER.out.vcf
    ch_versions = ch_versions.mix(FILTER.out.versions.first())


    // calculate pool statistics (fst, fisher, etc.)
    POOLSTATS(ch_filtered,ch_ref)
    ch_versions = ch_versions.mix(POOLSTATS.out.versions.first())

    // run post-processing steps
    POSTPROCESS(
        ch_vcf,
        POOLSTATS.out.fst,
        ch_ref,
        POOLSTATS.out.sync,
        POOLSTATS.out.split_sync,
        POOLSTATS.out.frequency
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  ''  + 'versions.yml',
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

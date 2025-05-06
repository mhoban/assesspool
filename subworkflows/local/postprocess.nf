include { GRENEDALF_FREQUENCY      } from '../../modules/local/grenedalf/frequency.nf'

workflow POSTPROCESS {

    take:
    ch_fst
    ch_ref
    ch_sync
    ch_split_sync

    main:

    ch_versions = Channel.empty()


    ch_ref = ch_ref.map{ meta, fasta, fai -> [ meta, fai ]}


    // ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // // TODO nf-core: edit emitted channels

    versions = ch_versions                     // channel: [ versions.yml ]
}


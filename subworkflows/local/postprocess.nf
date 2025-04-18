include { GRENEDALF_FREQUENCY      } from '../../modules/local/grenedalf/frequency.nf'

workflow PROCESSFST {

    take:
    ch_fst
    ch_ref
    ch_sync
    ch_split_sync

    main:

    ch_versions = Channel.empty()


    ch_ref = ch_ref.map{ meta, fasta, fai -> [ meta, fai ]}

    GRENEDALF_FREQUENCY(ch_sync,[[],[]],ch_ref)

    // ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // // TODO nf-core: edit emitted channels
    // bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    // bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    // csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}


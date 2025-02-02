include { PREPSYNC                              } from '../../modules/local/prepsync/main'
include { POPOOLATION2_FST                      } from '../../modules/local/popoolation2/fst.nf'
include { POPOOLATION2_FISHERTEST               } from '../../modules/local/popoolation2/fishertest.nf'
include { GRENEDALF_FST                         } from '../../modules/local/grenedalf/fst.nf'
include { POOLFSTAT_FST                         } from '../../modules/local/poolfstat/fst/main'

workflow FST {
    take:
    ch_vcf // channel: [ val(meta), path(vcf), path(index) ]
    ch_ref // channel:  path(ref)

    main:

    ch_versions = Channel.empty()

    ch_vcf
        .map{ meta, vcf, index -> [ meta.id, meta, vcf, index ] }
        .join( ch_ref.map{ meta, ref -> [ meta.id, ref ] } )
        .map{ id, meta, vcf, index ,ref -> [ meta, vcf, index, ref ] }
        .set{ ch_prep }

    // prepare sync file(s)
    PREPSYNC(ch_prep)

    // create split sync file channel
    PREPSYNC.out.split_sync
        .transpose()
        .map { meta, sync ->
            def poolsize = meta.pools.findAll { k, v -> sync =~ /\b${k}\b/}.sort { it.key }
            [ [ id: meta.id, pools: poolsize ], sync, poolsize ]
        }
        .set { ch_sync }

    // run popoolation fst if requested
    if (params.popoolation2) {
        POPOOLATION2_FST( ch_sync )
        ch_versions = ch_versions.mix(POPOOLATION2_FST.out.versions.first())

        if (params.fisher_test) {
            POPOOLATION2_FISHERTEST( ch_sync )
        }
    }

    // run grenedalf fst if requested
    if (params.grenedalf) {
        PREPSYNC.out.combined_sync
            .collectFile { meta, sync ->
                [ "${meta.id}.ps", meta.pools.collect{ k, v -> "${k}\t${v}"}.join("\n") ]
            }
            .map { f -> [ f.baseName, f ] }
            .set { ch_poolsizes }
        PREPSYNC.out.combined_sync
            .map { meta,sync -> [ meta.id, meta, sync ] }
            .join(ch_poolsizes)
            .map { id, meta, sync, poolsize -> [ meta, sync, poolsize ] }
            .set{ ch_gd }

        GRENEDALF_FST( ch_gd )
        ch_versions = ch_versions.mix(GRENEDALF_FST.out.versions.first())
    }

    // run poolfstat fst if requested
    if (params.poolfstat) {
        POOLFSTAT_FST( ch_sync )
        ch_versions = ch_versions.mix(POOLFSTAT_FST.out.versions.first())
    }


    emit:

    versions = ch_versions                     // channel: [ versions.yml ]
}


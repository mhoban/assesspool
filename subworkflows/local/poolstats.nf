include { POPOOLATION2_FST        } from '../../modules/local/popoolation2/fst.nf'
include { POPOOLATION2_FISHERTEST } from '../../modules/local/popoolation2/fishertest.nf'
include { GRENEDALF_FST           } from '../../modules/local/grenedalf/fst.nf'
include { GRENEDALF_SYNC          } from '../../modules/local/grenedalf/sync.nf'
include { POOLFSTAT_FST           } from '../../modules/local/poolfstat/fst/main'
include { GAWK as SPLIT_SYNC      } from '../../modules/nf-core/gawk/main'

/* function to produce a list of all possible combinations of m elements from a list */
def combn(list,m) {
    def n = list.size()
    m == 0 ?
        [[]] :
        (0..(n-m)).inject([]) { newlist, k ->
            def sublist = (k+1 == n) ? [] : list[(k+1)..<n]
            newlist += combn(sublist,m-1).collect { [list[k]] + it }
        }
}


workflow FST {
    take:
    ch_vcf // channel: [ val(meta), path(vcf), path(index) ]
    ch_ref // channel:  path(ref)

    main:

    ch_versions = Channel.empty()

    ch_vcf
        .map{ meta, vcf, index -> [ meta.id, meta, vcf, index ] }
        .join( ch_ref.map{ meta, ref, fai -> [ meta.id, ref, fai ] } )
        .map{ id, meta, vcf, index, ref, fai -> [ meta, vcf, index, ref, fai ] }
        .set{ ch_prep }

    // create sync file with grenedalf
    GRENEDALF_SYNC(
        ch_prep.map{ meta, vcf, index, ref, fai -> [ meta, vcf, index ] },
        ch_prep.map{ meta, vcf, index, ref, fai -> [ meta, [] ] },
        ch_prep.map{ meta, vcf, index, ref, fai -> [ meta, fai ] }
    )
    ch_versions = ch_versions.mix(GRENEDALF_SYNC.out.versions.first())


    // prepare and run either processes needing split sync files
    if (params.popoolation2 || params.poolfstat) {
        // split sync file into pairwise combinations of pools
        GRENEDALF_SYNC.out.sync
            .map{ meta, sync -> [ meta, sync, combn(meta.pools.keySet() as String[],2) ] }
            .transpose()
            .map{ meta, sync, pair ->
                [ [ id: meta.id, pools: meta.pools.subMap(pair).sort{ it.key } ], sync ]
            }
            .set{ ch_split_sync }

        SPLIT_SYNC( ch_split_sync, [] ).output
            .map{ meta, sync -> [ meta, sync, meta.pools ] }
            .set{ ch_split_sync }

        // run popoolation fst if requested
        if (params.popoolation2) {
            POPOOLATION2_FST( ch_split_sync )
            ch_versions = ch_versions.mix(POPOOLATION2_FST.out.versions.first())

            if (params.fisher_test) {
                POPOOLATION2_FISHERTEST( ch_split_sync )
            }
        }

        // run poolfstat fst if requested
        if (params.poolfstat) {
            POOLFSTAT_FST( ch_split_sync )
            ch_versions = ch_versions.mix(POOLFSTAT_FST.out.versions.first())
        }
    }

    // run grenedalf fst if requested
    if (params.grenedalf) {
        GRENEDALF_SYNC.out.sync
            .collectFile { meta, sync ->
                [ "${meta.id}.ps", meta.pools.collect{ k, v -> "${k}\t${v}"}.join("\n") ]
            }
            .map { f -> [ f.baseName, f ] }
            .set { ch_poolsizes }
        GRENEDALF_SYNC.out.sync
            .map { meta, sync -> [ meta.id, meta, sync ] }
            .join( ch_poolsizes )
            .join( ch_ref.map{ meta, ref, index -> [ meta.id, index ] } )
            .map { id, meta, sync, poolsize, index  -> [ meta, sync, poolsize, index ] }
            .set{ ch_gd }


        GRENEDALF_FST(
            ch_gd.map{ meta, sync, poolsize, index -> [ meta, sync, poolsize ] },
            ch_gd.map{ meta, sync, poolsize, index -> [ meta, [] ] },
            ch_gd.map{ meta, sync, poolsize, index -> [ meta, index ] }
        )

        // GRENEDALF_FST.out.fst.view()
        ch_versions = ch_versions.mix(GRENEDALF_FST.out.versions.first())
    }

    emit:

    versions = ch_versions                     // channel: [ versions.yml ]
}


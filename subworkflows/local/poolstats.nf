include { POPOOLATION2_FST        } from '../../modules/local/popoolation2/fst.nf'
include { POPOOLATION2_FISHERTEST } from '../../modules/local/popoolation2/fishertest.nf'
include { GRENEDALF_FST           } from '../../modules/local/grenedalf/fst.nf'
include { GRENEDALF_SYNC          } from '../../modules/local/grenedalf/sync.nf'
include { GRENEDALF_FREQUENCY     } from '../../modules/local/grenedalf/frequency.nf'
include { POOLFSTAT_FST           } from '../../modules/local/poolfstat/fst/main'
include { FISHERTEST              } from '../../modules/local/fishertest/main'
include { JOINFREQ                } from '../../modules/local/joinfreq/main'
include { GAWK as SPLIT_SYNC      } from '../../modules/nf-core/gawk/main'
include { GAWK as REHEADER_FST    } from '../../modules/nf-core/gawk/main'
include { GAWK as REHEADER_FISHER } from '../../modules/nf-core/gawk/main'
include { combn                   } from "../../modules/local/popoolation2/combn.nf"

workflow POOLSTATS {
    take:
    ch_input // channel: [ val(meta), path(vcf), path(index) ]
    ch_ref // channel:  path(ref)

    main:

    ch_split_sync = Channel.empty() // channel to receive split sync files
    ch_versions = Channel.empty()   // channel to receive tool versions
    ch_fst = Channel.empty()        // channel to receive fst output
    ch_fisher = Channel.empty()     // channel to receive fisher output

    // prepare vcf channel for input to grenedalf sync
    // (join fasta reference index)
    ch_input
        .map{ meta, file, index -> [ meta.id, meta, file, index ] }
        .join( ch_ref.map{ meta, ref, fai -> [ meta.id, ref, fai ] } )
        .map{ id, meta, file, index, ref, fai -> [ meta, file, index, ref, fai ] }
        .set{ ch_prep }

    ch_prep
        .map{ meta, file, index, ref, fai -> [ meta, meta.rename ? meta.pool_map : [:] ] }
        // .filter { it[0].rename }
        .collectFile(sort: true) { meta, pool_map -> [ "${meta.id}.map", pool_map ? pool_map.collect{k,v -> "${k}\t${v}"}.join("\n") : "" ] }
        .map{ [ it.baseName, it ]}
        .set{ ch_samplemap }

    ch_prep
        .map{ [ it[0].id ] + it  }
        .join( ch_samplemap )
        .map{ id, meta, file, index, ref, fai, pool_map ->
            [meta, file, index, ref, fai, pool_map.size() ? pool_map : [] ]
        }
        .set{ ch_prep }


    // create/modify sync files using grenedalf
    GRENEDALF_SYNC(
        ch_prep.map{ meta, file, index, ref, fai, map -> [ meta, meta.format == 'vcf' ? file : [], index ] },
        ch_prep.map{ meta, file, index, ref, fai, map -> [ meta, meta.format == 'sync' ? file : [] ] },
        ch_prep.map{ meta, file, index, ref, fai, map -> [ meta, [] ] },
        ch_prep.map{ meta, file, index, ref, fai, map -> [ meta, fai ] },
        ch_prep.map{ meta, file, index, ref, fai, map -> [ meta, map ] }
    )
    ch_versions = ch_versions.mix(GRENEDALF_SYNC.out.versions.first())

    // prepare and split input channel into pairwise comparisons
    if (params.popoolation2 || params.poolfstat || params.fisher_test_popoolation) {
        GRENEDALF_SYNC.out.sync
            .map{ meta, sync -> [ meta, sync, combn(meta.pools.keySet() as String[],2) ] }
            .transpose()
            .map{ meta, sync, pair ->
                // println meta.pools.subMap(pair).sort{ it.key }.values().toArray().getClass()
                [
                    meta + [
                        pools: meta.pools.subMap(pair).sort{ it.key },
                        pool_map: meta.pool_map.findAll{ it.value as String in pair as String[] }.sort{ it.key }
                    ],
                    sync
                ]
            }
            .set{ ch_split_sync }

        // pairwise split sync file
        SPLIT_SYNC( ch_split_sync, [] ).output
            .map{ meta, sync -> [ meta, meta.pools, sync ] }
            .set{ ch_split_sync }
        ch_versions = ch_versions.mix(SPLIT_SYNC.out.versions.first())
    }

    // prepare input for frequency calculation
    GRENEDALF_SYNC.out.sync
        .map { meta, sync -> [ meta.id, meta, sync ] }
        .join( ch_ref.map{ meta, ref, index -> [ meta.id, index ] } )
        .map { id, meta, sync, index  -> [ meta, sync, index ] }
        .set{ ch_freq }

    // calculate snp frequencies using grenedalf
    ch_freq = GRENEDALF_FREQUENCY(
        ch_freq.map{ meta, sync, index -> [ meta, sync ] },
        ch_freq.map{ meta, sync, index -> [ meta, [] ] },
        ch_freq.map{ meta, sync, index -> [ meta, index ] }
    ).freq
    ch_versions = ch_versions.mix(GRENEDALF_FREQUENCY.out.versions.first())

    // run popoolation fst if requested
    if (params.popoolation2) {
        POPOOLATION2_FST( ch_split_sync )
        REHEADER_FST( POPOOLATION2_FST.out.fst, [] )
        ch_versions = ch_versions.mix(POPOOLATION2_FST.out.versions.first())
        ch_versions = ch_versions.mix(REHEADER_FST.out.versions.first())

        ch_fst = ch_fst.mix( REHEADER_FST.out.output.map { meta, fst -> [ meta, fst, "popoolation" ] } )
    }

    // run poolfstat fst if requested
    if (params.poolfstat) {
        POOLFSTAT_FST( ch_split_sync )
        ch_versions = ch_versions.mix(POOLFSTAT_FST.out.versions.first())

        // mix fst results
        ch_fst = ch_fst.mix( POOLFSTAT_FST.out.fst.map { meta, fst -> [ meta, fst, "poolfstat" ] }  )
    }

    // run grenedalf fst if requested
    if (params.grenedalf) {
        // create pool size map file
        GRENEDALF_SYNC.out.sync
            .collectFile { meta, sync ->
                [ "${meta.id}.ps", meta.pools.collect{ k, v -> "${k}\t${v}"}.join("\n") ]
            }
            .map { f -> [ f.baseName, f ] }
            .set { ch_poolsizes }

        // create grenedalf input channel
        GRENEDALF_SYNC.out.sync
            .map { meta, sync -> [ meta.id, meta, sync ] }
            .join( ch_poolsizes )
            .join( ch_ref.map{ meta, ref, index -> [ meta.id, index ] } )
            .map { id, meta, sync, poolsize, index  -> [ meta, sync, poolsize, index ] }
            .set{ ch_gd }

        // run grenedalf fst calculation
        GRENEDALF_FST(
            ch_gd.map{ meta, sync, poolsize, index -> [ meta, sync, poolsize ] },
            ch_gd.map{ meta, sync, poolsize, index -> [ meta, [] ] },
            ch_gd.map{ meta, sync, poolsize, index -> [ meta, index ] }
        )
        ch_versions = ch_versions.mix(GRENEDALF_FST.out.versions.first())

        // mix fst results
        ch_fst = ch_fst.mix( GRENEDALF_FST.out.fst.map { meta, fst -> [ meta, fst, "grenedalf" ] }  )
    }

    // prepare input channel to join fst results with snp frequencies
    ch_fst
        .map { meta, fst, method -> [ meta.id, meta, fst, method ] }
        .combine( ch_freq.map { meta, freq -> [ meta.id, freq ] }, by: 0 )
        .map { id, meta, fst, method, freq -> [ meta, freq, fst, method ] }
        .set { ch_join }

    // join pairwise fst results to snp frequencies and concatenate into one big file per input
    ch_join = JOINFREQ( ch_join ).fst_freq
        .collectFile( keepHeader: true, sort: true, storeDir: 'output/fst', cache: true ){ meta, joined -> [ "${meta.id}.fst", joined ] }
        .map{ [it.baseName, it ] }

    // rejoin meta tags
    ch_join = ch_input
        .map{ it[0] }
        .unique()
        .map{ [ it.id, it ] }
        .join( ch_join )
        .map{ id, meta, f -> [ meta, f ] }

    ch_versions = ch_versions.mix(JOINFREQ.out.versions.first())

    // run fisher tests if requested
    if (params.fisher_test) {
        // pull out population pairs from frequency table
        ch_freq
            .map{ meta, freq -> [ meta, freq, combn(meta.pools.keySet() as String[],2) ] }
            .transpose()
            .map{ meta, freq, pair ->
                [
                    meta + [
                        pools: meta.pools.subMap(pair).sort{ it.key },
                        pool_map: meta.pool_map.findAll{ it.value as String in pair as String[] }.sort{ it.key }
                    ],
                    pair.sort(),
                    freq
                ]
            }
            .set{ ch_split_freq }

        // run pairwise fisher tests
        FISHERTEST( ch_split_freq ).fisher

        // mix fisher results & versions
        ch_fisher = ch_fisher.mix( FISHERTEST.out.fisher )
        ch_versions = ch_versions.mix(FISHERTEST.out.versions.first())

    }

    // run fisher tests with popoolation2
    if (params.fisher_test_popoolation) {
        POPOOLATION2_FISHERTEST( ch_split_sync )
        REHEADER_FISHER( POPOOLATION2_FISHERTEST.out.fisher, [] )

        // mix fisher results & versions
        // ch_fisher = ch_fisher.mix( POPOOLATION2_FISHERTEST.out.fisher )
        ch_fisher = ch_fisher.mix( REHEADER_FISHER.out.output )
        ch_versions = ch_versions.mix(POPOOLATION2_FISHERTEST.out.versions.first())
        ch_versions = ch_versions.mix(REHEADER_FISHER.out.versions.first())
    }


    emit:
    fst = ch_join
    all_fst = ch_fst
    fisher = ch_fisher
    sync = GRENEDALF_SYNC.out.sync
    frequency = ch_freq
    split_sync = ch_split_sync
    versions = ch_versions                     // channel: [ versions.yml ]
}


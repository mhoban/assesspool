
    take:
    ch_vcf // channel: [ val(meta), path(vcf) ]
    ch_ref // channel:  path(ref)

    main:

    ch_versions = Channel.empty()

    // prepare sync file(s)
    PREPSYNC(ch_vcf,ch_ref)

    // create split sync file channel
    PREPSYNC.out.split_sync
        .transpose()
        .map { meta, sync ->
            [ [ id: meta.id, populations: meta.populations.findAll { k, v -> sync =~ /\b${k}\b/}.sort { it.key } ], sync ]
        }
        .set { ch_split_sync }

    // run popoolation fst if requested
    if (params.popoolation2) {
        POPOOLATION2_FST( ch_split_sync )
        ch_versions = ch_versions.mix(POPOOLATION2_FST.out.versions.first())
    }

    // run grenedalf fst if requested
    if (params.grenedalf) {

        PREPSYNC.out.combined_sync
            .map { meta, sync -> meta.populations.collect { k, v -> "${k}\t${v}\n" } }
            .flatten()
            .collectFile(name: "pool_sizes.txt", cache: false)
            .set { ch_poolsizes }

        GRENEDALF_FST( PREPSYNC.out.combined_sync, ch_poolsizes )
        ch_versions = ch_versions.mix(GRENEDALF_FST.out.versions.first())
    }

    // run poolfstat fst if requested
    if (params.poolfstat) {
        POOLFSTAT_FST( ch_split_sync )
        ch_versions = ch_versions.mix(POOLFSTAT_FST.out.versions.first())
    }


    emit:

    versions = ch_versions                     // channel: [ versions.yml ]
}


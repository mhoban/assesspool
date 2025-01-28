include { BCFTOOLS_QUERY as VCF_SAMPLES   } from '../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_REHEADER as VCF_RENAME } from '../../modules/nf-core/bcftools/reheader/main'

workflow PREPROCESS {

    take:
    ch_samplesheet

    main:

    ch_versions = Channel.empty()

    // save ref genome channel
    ch_samplesheet
        .map { meta, vcf, ref -> ref }
        .set { ch_ref }

    // save vcf channel
    ch_samplesheet
        .map { meta, vcf, ref -> [ meta, vcf ] }
        .set { ch_vcf }


    // Get VCF sample names from header
    VCF_SAMPLES(ch_vcf.combine(Channel.fromPath('fake.tbi')),[],[],[]).output
        .splitText()
        .map { meta, pool -> [ meta.id, pool.trim() ] }
        .groupTuple()
        .set { ch_samplenames }
    ch_versions = ch_versions.mix(VCF_SAMPLES.out.versions.first())


    // Extract VCF sample names and smash them into a pool map
    ch_vcf
        .map{ meta, vcf -> [ meta.id, meta, vcf ] }
        .join(ch_samplenames)
        .map { id, meta, vcf, pools ->
            def pp = meta.pools ?: pools
            def ps = meta.pool_sizes.size() > 1 ? meta.pool_sizes : (1..pools.size()).collect{ meta.pool_sizes[0] }
            [ [ id: meta.id, rename: meta.rename, pools: [pp,ps].transpose().collectEntries(), pool_map: [pools,pp].transpose().collectEntries() ], vcf ]
        }
        .set{ ch_vcf }


    // Rename any VCF headers we were asked to
    // ---------------------------------------
    // First branch out which ones we need to rename
    ch_vcf
        .branch { meta, vcf ->
            rename: meta.rename
            leave: !meta.rename
        }
        .set{ ch_vcf }

    // Create a sample name map file for bcftools reheader
    ch_vcf.rename
        // mapfiles named for meta.id
        .collectFile { meta, vcf ->
            [ "${meta.id}.map", meta.pool_map.collect{ k, v -> "${k} ${v}"}.join("\n") ]
        }
        // use mapfile basename for meta.id to join below
        .map { f -> [ f.baseName, f ] }
        .set{ ch_remap }

    // Join sample map back to VCF channel
    ch_vcf.rename
        .map { meta, vcf -> [ meta.id, meta, vcf ] }
        .join(ch_remap)
        .map { id, meta, vcf, mapfile -> [ meta, vcf, [], mapfile ] }
        .set{ ch_remap }

    // Do the renaming we need to do
    VCF_RENAME(ch_remap,ch_remap.map{ [it[0],[]]}).vcf
        .set { ch_remap }
    // mix in versions.yml
    ch_versions = ch_versions.mix(VCF_RENAME.out.versions.first())
    // ---------------------------------------

    // Mix back in any renamed VCFs
    ch_vcf.leave
        .mix(ch_remap)
        .set{ ch_vcf }


    emit:
    vcf = ch_vcf
    ref = ch_ref

    versions = ch_versions                     // channel: [ versions.yml ]
}


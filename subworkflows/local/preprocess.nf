include { BCFTOOLS_QUERY as VCF_SAMPLES     } from '../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_REHEADER as VCF_RENAME   } from '../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_VIEW as COMPRESS_VCF     } from '../../modules/nf-core/bcftools/view/main'
include { TABIX_TABIX as INDEX_VCF          } from '../../modules/nf-core/tabix/tabix/main'
include { SAMTOOLS_FAIDX as INDEX_REFERENCE } from '../../modules/nf-core/samtools/faidx/main'

workflow PREPROCESS {

    take:
    ch_samplesheet

    main:

    ch_versions = Channel.empty()

    // save ref genome channel
    ch_samplesheet
        .map { meta, vcf, index, ref -> [ meta, ref ] }
        .set { ch_ref }
    // generate fasta ref index from reference fasta
    INDEX_REFERENCE( ch_ref, ch_ref.map{ meta, ref -> [ meta, [] ] } )

    ch_ref
        .map{ meta, fasta -> [ meta.id, meta, fasta ] }
        .join( INDEX_REFERENCE.out.fai.map{ meta, fai -> [ meta.id, fai ] } )
        .map{ id, meta, fasta, fai -> [ meta, fasta, fai ] }
        .set{ ch_ref }

    // save vcf channel
    ch_samplesheet
        .map { meta, vcf, index, ref -> [ meta, vcf, index ] }
        .set { ch_vcf }


    // Get VCF sample names from header regardless of whether we're renaming them
    VCF_SAMPLES(ch_vcf,[],[],[]).output
        .splitText()
        .map { meta, pool -> [ meta.id, pool.trim() ] }
        .groupTuple()
        .set { ch_samplenames }
    ch_versions = ch_versions.mix(VCF_SAMPLES.out.versions.first())


    // Extract VCF sample names and smash them into a pool map
    ch_vcf
        .map{ meta, vcf, index -> [ meta.id, meta, vcf, index ] }
        .join( ch_samplenames )
        .map { id, meta, vcf, index, pools ->
            def pp = meta.pools ?: pools
            def ps = meta.pool_sizes.size() > 1 ? meta.pool_sizes : (1..pools.size()).collect{ meta.pool_sizes[0] }
            [ [ id: meta.id, rename: meta.rename, pools: [pp,ps].transpose().collectEntries(), pool_map: [pools,pp].transpose().collectEntries() ], vcf, index ]
        }
        .set{ ch_vcf }

    // First, branch out which ones we need to rename, zip or index
    ch_vcf
        .branch { meta, vcf, index ->
            rename: meta.rename
            zip: !meta.rename && vcf.extension != "gz"
            index: !meta.rename && vcf.extension == "gz" && !index
            keep_as_is: !meta.rename && vcf.extension == "gz" && index
        }
        .set{ ch_vcf }

    // rename VCF samples -----------------------------------
    // Create a sample name map file for bcftools reheader
    ch_vcf.rename
        // mapfiles named for meta.id
        .collectFile { meta, vcf, index ->
            [ "${meta.id}.map", meta.pool_map.collect{ k, v -> "${k} ${v}"}.join("\n") ]
        }
        // use mapfile basename for meta.id to join below
        .map { f -> [ f.baseName, f ] }
        .set{ ch_remap }

    // Join sample map back to VCF channel
    ch_vcf.rename
        .map { meta, vcf, index -> [ meta.id, meta, vcf ] }
        .join( ch_remap )
        .map { id, meta, vcf, mapfile -> [ meta, vcf, [], mapfile ] }
        .set{ ch_remap }

    // Rename samples in the VCF header
    VCF_RENAME( ch_remap,ch_remap.map{ [it[0],[]]} ).vcf
        .map{ meta, vcf -> [meta.id, meta, vcf] }
        .join( VCF_RENAME.out.index.map{ meta, index -> [ meta.id, index ] } )
        .map{ id, meta, vcf, index -> [ meta, vcf, index ] }
        .set { ch_renamed }
    ch_versions = ch_versions.mix(VCF_RENAME.out.versions.first())
    // rename VCF samples -----------------------------------

    // zip unzipped VCF -------------------------------------
    COMPRESS_VCF( ch_vcf.zip,[],[],[] )
    COMPRESS_VCF.out.vcf
        .map{ meta,vcf -> [ meta.id, meta, vcf ] }
        .join( COMPRESS_VCF.out.tbi.map{ meta,index -> [ meta.id, index ] } )
        .map{ id, meta, vcf, index -> [meta,vcf,index] }
        .set{ ch_zipped }
    ch_versions = ch_versions.mix(COMPRESS_VCF.out.versions.first())
    // zip unzipped VCF -------------------------------------

    // index unindexed VCF ----------------------------------
    INDEX_VCF( ch_vcf.index.map{ meta, vcf, index -> [ meta,vcf ] } )
    ch_vcf.index.map{ meta, vcf, index -> [ meta.id, meta, vcf ] }
        .join(INDEX_VCF.out.tbi.map{ meta, index -> [ meta.id, index ] } )
        .map{ id, meta, vcf, index -> [meta, vcf, index] }
        .set{ ch_indexed }
    ch_versions = ch_versions.mix(INDEX_VCF.out.versions.first())
    // index unindexed VCF ----------------------------------


    // // Mix back in any renamed tbis
    ch_vcf.keep_as_is
        .mix( ch_renamed )
        .mix( ch_zipped )
        .mix( ch_indexed )
        .set{ ch_vcf }


    emit:
    vcf = ch_vcf
    ref = ch_ref

    versions = ch_versions                     // channel: [ versions.yml ]
}


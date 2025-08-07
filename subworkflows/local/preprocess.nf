include { BCFTOOLS_QUERY as VCF_SAMPLES     } from '../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_VIEW as COMPRESS_VCF     } from '../../modules/nf-core/bcftools/view/main'
include { TABIX_TABIX as INDEX_VCF          } from '../../modules/nf-core/tabix/tabix/main'
include { SAMTOOLS_FAIDX as INDEX_REFERENCE } from '../../modules/nf-core/samtools/faidx/main'

// read the first line from a gzip'd file
def gz_head(Path file) {
    new java.io.BufferedReader(
        new java.io.InputStreamReader(
            new java.util.zip.GZIPInputStream(new java.io.FileInputStream(file.toFile()))
        )
    ).withCloseable { reader ->
        return reader.readLine()
    }
}

workflow PREPROCESS {

    take:
    ch_samplesheet

    main:

    ch_versions = Channel.empty()

    // save ref genome channel
    ch_samplesheet
        .map { meta, input, index, ref -> [ meta, ref ] }
        .set { ch_ref }
    // generate fasta ref index from reference fasta
    INDEX_REFERENCE( ch_ref, ch_ref.map{ meta, ref -> [ meta, [] ] } )

    ch_ref
        .map{ meta, fasta -> [ meta.id, meta, fasta ] }
        .join( INDEX_REFERENCE.out.fai.map{ meta, fai -> [ meta.id, fai ] } )
        .map{ id, meta, fasta, fai -> [ meta, fasta, fai ] }
        .set{ ch_ref }

    // save input channel
    ch_samplesheet
        .map { meta, input, index, ref -> [ meta, input, index ] }
        .set { ch_input }


    // branch inputs into sync and/or vcf
    ch_input
        .branch { meta, input, index ->
            vcf: input.toString() =~ /(?i)\.vcf(\.gz)?$/
            sync: input.toString() =~ /(?i)\.sync(\.gz)?$/
        }
        .set{ ch_input }

    // preprocess VCF input(s)
    // Get VCF sample names from header regardless of whether we're renaming them
    VCF_SAMPLES(ch_input.vcf,[],[],[]).output
        .splitText()
        .map { meta, pool -> [ meta.id, pool.trim() ] }
        .groupTuple()
        .set { ch_samplenames }

    ch_versions = ch_versions.mix(VCF_SAMPLES.out.versions.first())

    // Extract VCF sample names and smash them into a pool map
    ch_input.vcf
        .map{ meta, vcf, index -> [ meta.id, meta, vcf, index ] }
        .join( ch_samplenames )
        .map { id, meta, vcf, index, pools ->
            def pp = meta.pools ?: pools
            def ps = meta.pool_sizes.size() > 1 ? meta.pool_sizes : (1..pools.size()).collect{ meta.pool_sizes[0] }
            [ [ id: meta.id, format: 'vcf', rename: meta.rename, pools: [pp,ps].transpose().collectEntries(), pool_map: [pools,pp].transpose().collectEntries() ], vcf, index ]
        }
        .set{ ch_vcf }

    // preprocess sync input(s)
    // Extract/get sync sample names the same way, more or less
	ch_input.sync
        .map{ meta, sync, index ->
            def pools = (sync.extension.toLowerCase() == 'gz' ? gz_head(sync) : sync.withReader { it.readLine() }).split(/\t/)
            def sync_hdr = !!(pools[0] =~ /^#chr/)
            pools = sync_hdr ? pools[3..-1] : (1..(pools.size()-3)).collect{ "${sync.baseName}.${it}"}

            def pp = meta.pools ?: pools
            assert pp.size() == pools.size()
            def ps = meta.pool_sizes.size() > 1 ? meta.pool_sizes : (1..pools.size()).collect{ meta.pool_sizes[0] }
            [
                meta + [
                    format: 'sync',
                    sync_hdr: sync_hdr,
                    pools: [pp,ps].transpose().collectEntries{k,v -> [ k.toString(), v ]},
                    pool_map: [pools,pp].transpose().collectEntries{k,v -> [ k.toString(), v.toString() ]}
                ],
                sync,
                []
            ]
        }
        .set { ch_sync }

    // Figure out how to process vcf samples
    ch_vcf
        .branch { meta, vcf, index ->
            // rename: meta.rename
            zip: /*!meta.rename &&*/ vcf.extension != "gz"
            index: /*!meta.rename &&*/ vcf.extension == "gz" && !index
            keep_as_is: /*!meta.rename &&*/ vcf.extension == "gz" && index
        }
        .set{ ch_vcf }


    // zip unzipped VCF -------------------------------------
    COMPRESS_VCF( ch_vcf.zip,[],[],[] )
    COMPRESS_VCF.out.vcf
        .map{ meta,vcf -> [ meta.id, meta, vcf ] }
        .join( COMPRESS_VCF.out.tbi.map{ meta,index -> [ meta.id, index ] } )
        .map{ id, meta, vcf, index -> [meta,vcf,index] }
        .set{ ch_zipped }
    ch_versions = ch_versions.mix(COMPRESS_VCF.out.versions.first())

    // index unindexed VCF ----------------------------------
    INDEX_VCF( ch_vcf.index.map{ meta, vcf, index -> [ meta,vcf ] } )
    ch_vcf.index.map{ meta, vcf, index -> [ meta.id, meta, vcf ] }
        .join(INDEX_VCF.out.tbi.map{ meta, index -> [ meta.id, index ] } )
        .map{ id, meta, vcf, index -> [meta, vcf, index] }
        .set{ ch_indexed }
    ch_versions = ch_versions.mix(INDEX_VCF.out.versions.first())


    // Mix back in any renamed bits
    ch_vcf.keep_as_is
        .mix( ch_zipped )
        .mix( ch_indexed )
        .mix( ch_sync )
        .set{ ch_input }

    emit:
    output = ch_input
    ref = ch_ref

    versions = ch_versions                     // channel: [ versions.yml ]
}


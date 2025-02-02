include { VCFTOOLS as VCFTOOLS_FILTER                    } from '../../modules/nf-core/vcftools/main'
include { VCFTOOLS as VCFTOOLS_THIN                      } from '../../modules/nf-core/vcftools/main'
include { BCFTOOLS_FILTER                                } from '../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_VIEW as BCFTOOLS_COMPRESS_INDEX_FILTERED } from '../../modules/nf-core/bcftools/view/main'

workflow FILTER {

    take:
    ch_vcf // channel: [ val(meta), [ vcf ] ]

    main:

    ch_versions = Channel.empty()


    // decide whether to run bcftools filter
    def bcff = params.match_allele_count || ( params.filter && ( params.min_mapping_quality || params.min_mapping_quality_ref ||
                params.min_mapping_ratio || params.max_mapping_ratio || params.read_balance_left || params.read_balance_right ||
                params.quality_depth_ratio || params.mispaired_reads || params.min_pools ||
                params.min_depth || params.max_allele_length || params.min_quality ||
                params.variant_type || params.min_alternate_observations ) )

    if (bcff) {
        BCFTOOLS_FILTER( ch_vcf )
        BCFTOOLS_FILTER.out.vcf
            .map{ meta, vcf -> [ meta.id, meta, vcf ] }
            .join( BCFTOOLS_FILTER.out.tbi.map{ meta, index -> [ meta.id, index ] } )
            .map{ id, meta, vcf, index -> [ meta, vcf, index ] }
            .set{ ch_vcf2 }
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())
    } else {
        ch_vcf2 = ch_vcf
    }

    // decide whether to run vcftools filter
    def vcft = params.filter && ( params.max_missing || params.min_minor_allele_count ||params.min_mean_depth ||
            params.max_mean_depth || params.hwe_cutoff )

    if (vcft) {
        VCFTOOLS_FILTER ( ch_vcf2.map{ meta, vcf, index -> [ meta, vcf ] }, [], [] )
        ch_versions = ch_versions.mix(VCFTOOLS_FILTER.out.versions.first())
        // TODO: re-compress and index output
        VCFTOOLS_FILTER.out.vcf
            .map{ meta, vcf -> [ meta, vcf, [] ] }
            .set{ ch_vcf3 }
    } else {
        ch_vcf3 = ch_vcf2
    }

    if (params.thin_snps) {
        VCFTOOLS_THIN( ch_vcf3.map{ meta, vcf, index -> [ meta, vcf ] }, [], [] )
        ch_versions = ch_versions.mix(VCFTOOLS_THIN.out.versions.first())

        VCFTOOLS_THIN.out.vcf
            .map{ meta, vcf -> [ meta, vcf, [] ] }
            .set{ ch_vcf4 }
    } else {
        ch_vcf4 = ch_vcf3
    }

    if (vcft || params.thin_snps) {
        // re-index and compress if vcftools was run
        BCFTOOLS_COMPRESS_INDEX_FILTERED( ch_vcf4, [], [], [] ).vcf
            .map{ meta,vcf -> [ meta.id, meta, vcf ] }
            .join( BCFTOOLS_COMPRESS_INDEX_FILTERED.out.tbi.map{ meta,index -> [ meta.id, index ] } )
            .map{ id, meta, vcf, index -> [ meta, vcf, index ] }
            .set{ ch_vcf_final }
        ch_versions = ch_versions.mix(BCFTOOLS_COMPRESS_INDEX_FILTERED.out.versions.first())
    } else {
        ch_vcf_final = ch_vcf4
    }

    // TODO: visualize filter results? (rmd section 'filter_visualization')
    emit:
    vcf      = ch_vcf_final
    versions = ch_versions                     // channel: [ versions.yml ]
}


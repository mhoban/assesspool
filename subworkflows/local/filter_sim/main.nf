include { VCFTOOLS as MAX_MEAN_DEPTH                 } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as MIN_MEAN_DEPTH                 } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as HWE_CUTOFF                     } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as MINOR_ALLELE_COUNT             } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as MAX_MISSING                    } from '../../../modules/nf-core/vcftools/main'
include { BCFTOOLS_FILTER as MAX_ALLELE_LENGTH       } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as QUALITY_DEPTH_RATIO     } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as MIN_QUALITY             } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as VARIANT_TYPE            } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as MISPAIRED_READS         } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as ALTERNATE_OBSERVATIONS  } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as MIN_MAPPING_QUALITY     } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as MAPPING_RATIO           } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as MIN_DEPTH               } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as MIN_POOLS               } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as READ_BALANCE            } from '../../../modules/nf-core/bcftools/filter/main'
include { VCFTOOLS as THIN                           } from '../../../modules/nf-core/vcftools/main'
include { RIPGREP as COUNT_SNPS                      } from '../../../modules/nf-core/ripgrep/main'

include { BCFTOOLS_VIEW as BCFTOOLS_COMPRESS_INDEX_FILTERED } from '../../../modules/nf-core/bcftools/view/main'

workflow FILTER_SIM {

    take:
    ch_vcf

    main:
    ch_versions = Channel.empty()

    ch_filter_summary = ch_vcf.map{ meta, vcf, index -> [ meta + [filter: 'before'], vcf, index ] }

    // run vcftools-based filters
    if (params.max_mean_depth) {
        MAX_MEAN_DEPTH( ch_vcf.map { meta, vcf, index -> [ meta, vcf ] }, [], [] )
        ch_filter_summary = ch_filter_summary.mix( MAX_MEAN_DEPTH.out.vcf.map{ meta, vcf -> [ meta + [filter: 'max_mean_depth'], vcf, [] ] } )
    }

    if (params.min_mean_depth) {
        MIN_MEAN_DEPTH( ch_vcf.map { meta, vcf, index -> [ meta, vcf ] }, [], [] )
        ch_filter_summary = ch_filter_summary.mix( MIN_MEAN_DEPTH.out.vcf.map{ meta, vcf -> [ meta + [filter: 'min_mean_depth'], vcf, [] ] } )
    }

    if (params.hwe_cutoff) {
        HWE_CUTOFF( ch_vcf.map { meta, vcf, index -> [ meta, vcf ] }, [], [] )
        ch_filter_summary = ch_filter_summary.mix( HWE_CUTOFF.out.vcf.map{ meta, vcf -> [ meta + [filter: 'hwe_cutoff'], vcf, [] ] } )
    }

    if (params.min_minor_allele_count) {
        MINOR_ALLELE_COUNT( ch_vcf.map { meta, vcf, index -> [ meta, vcf ] }, [], [] )
        ch_filter_summary = ch_filter_summary.mix( MINOR_ALLELE_COUNT.out.vcf.map{ meta, vcf -> [ meta + [filter: 'minor_allele_count'], vcf, [] ] } )
    }

    if (params.max_missing) {
        MAX_MISSING( ch_vcf.map { meta, vcf, index -> [ meta, vcf ] }, [], [] )
        ch_filter_summary = ch_filter_summary.mix( MAX_MISSING.out.vcf.map{ meta, vcf -> [ meta + [filter: 'max_missing'], vcf, [] ] } )
    }

    if (params.thin_snps) {
        THIN( ch_vcf.map { meta, vcf, index -> [ meta, vcf ] }, [], [] )
        ch_filter_summary = ch_filter_summary.mix( THIN.out.vcf.map{ meta, vcf -> [ meta + [filter: 'thin'], vcf, [] ] } )
    }

    // run bcftools-based filters
    if (params.max_allele_length) {
        MAX_ALLELE_LENGTH( ch_vcf )
        ch_filter_summary = ch_filter_summary.mix(
            MAX_ALLELE_LENGTH.out.vcf
                .map{ meta, vcf -> [ meta.id, meta, vcf ] }
                .join( MAX_ALLELE_LENGTH.out.tbi.map{ meta, index -> [ meta.id, index ] } )
                .map{ id, meta, vcf, index -> [ meta + [filter: 'max_allele_length'], vcf, index ] }
        )
    }

    if (params.quality_depth_ratio) {
        QUALITY_DEPTH_RATIO( ch_vcf )
        ch_filter_summary = ch_filter_summary.mix(
            QUALITY_DEPTH_RATIO.out.vcf
                .map{ meta, vcf -> [ meta.id, meta, vcf ] }
                .join( QUALITY_DEPTH_RATIO.out.tbi.map{ meta, index -> [ meta.id, index ] } )
                .map{ id, meta, vcf, index -> [ meta + [filter: 'quality_depth_ratio'], vcf, index ] }
        )
    }

    if (params.min_quality) {
        MIN_QUALITY( ch_vcf )
        ch_filter_summary = ch_filter_summary.mix(
            MIN_QUALITY.out.vcf
                .map{ meta, vcf -> [ meta.id, meta, vcf ] }
                .join( MIN_QUALITY.out.tbi.map{ meta, index -> [ meta.id, index ] } )
                .map{ id, meta, vcf, index -> [ meta + [filter: 'min_quality'], vcf, index ] }
        )
    }

    if (params.variant_type) {
        VARIANT_TYPE( ch_vcf )
        ch_filter_summary = ch_filter_summary.mix(
            VARIANT_TYPE.out.vcf
                .map{ meta, vcf -> [ meta.id, meta, vcf ] }
                .join( VARIANT_TYPE.out.tbi.map{ meta, index -> [ meta.id, index ] } )
                .map{ id, meta, vcf, index -> [ meta + [filter: 'variant_type'], vcf, index ] }
        )
    }

    if (params.mispaired_reads) {
        MISPAIRED_READS( ch_vcf )
        ch_filter_summary = ch_filter_summary.mix(
            MISPAIRED_READS.out.vcf
                .map{ meta, vcf -> [ meta.id, meta, vcf ] }
                .join( MISPAIRED_READS.out.tbi.map{ meta, index -> [ meta.id, index ] } )
                .map{ id, meta, vcf, index -> [ meta + [filter: 'mispaired_reads'], vcf, index ] }
        )
    }

    if (params.min_alternate_observations) {
        ALTERNATE_OBSERVATIONS( ch_vcf )
        ch_filter_summary = ch_filter_summary.mix(
            ALTERNATE_OBSERVATIONS.out.vcf
                .map{ meta, vcf -> [ meta.id, meta, vcf ] }
                .join( ALTERNATE_OBSERVATIONS.out.tbi.map{ meta, index -> [ meta.id, index ] } )
                .map{ id, meta, vcf, index -> [ meta + [filter: 'alternate_observations'], vcf, index ] }
        )
    }


    if (params.min_mapping_quality) {
        MIN_MAPPING_QUALITY( ch_vcf )
        ch_filter_summary = ch_filter_summary.mix(
            MIN_MAPPING_QUALITY.out.vcf
                .map{ meta, vcf -> [ meta.id, meta, vcf ] }
                .join( MIN_MAPPING_QUALITY.out.tbi.map{ meta, index -> [ meta.id, index ] } )
                .map{ id, meta, vcf, index -> [ meta + [filter: 'min_mapping_quality'], vcf, index ] }
        )
    }

    if (params.min_mapping_ratio || params.max_mapping_ratio) {
        MAPPING_RATIO( ch_vcf )
        ch_filter_summary = ch_filter_summary.mix(
            MAPPING_RATIO.out.vcf
                .map{ meta, vcf -> [ meta.id, meta, vcf ] }
                .join( MAPPING_RATIO.out.tbi.map{ meta, index -> [ meta.id, index ] } )
                .map{ id, meta, vcf, index -> [ meta + [filter: 'mapping_ratio'], vcf, index ] }
        )
    }

    if (params.min_depth) {
        MIN_DEPTH( ch_vcf )
        ch_filter_summary = ch_filter_summary.mix(
            MIN_DEPTH.out.vcf
                .map{ meta, vcf -> [ meta.id, meta, vcf ] }
                .join( MIN_DEPTH.out.tbi.map{ meta, index -> [ meta.id, index ] } )
                .map{ id, meta, vcf, index -> [ meta + [filter: 'min_depth'], vcf, index ] }
        )
    }

    if (params.min_pools) {
        MIN_POOLS( ch_vcf )
        ch_filter_summary = ch_filter_summary.mix(
            MIN_POOLS.out.vcf
                .map{ meta, vcf -> [ meta.id, meta, vcf ] }
                .join( MIN_POOLS.out.tbi.map{ meta, index -> [ meta.id, index ] } )
                .map{ id, meta, vcf, index -> [ meta + [filter: 'min_pools'], vcf, index ] }
        )
    }

    if (params.read_balance_left || params.read_balance_right) {
        READ_BALANCE( ch_vcf )
        ch_filter_summary = ch_filter_summary.mix(
            READ_BALANCE.out.vcf
                .map{ meta, vcf -> [ meta.id, meta, vcf ] }
                .join( READ_BALANCE.out.tbi.map{ meta, index -> [ meta.id, index ] } )
                .map{ id, meta, vcf, index -> [ meta + [filter: 'read_balance'], vcf, index ] }
        )
    }


    COUNT_SNPS( ch_filter_summary.map{ meta, vcf, index -> [ meta, vcf ] }, "'^#'", false )
    ch_snp_counts = COUNT_SNPS.out.txt
    ch_snp_counts.map{ meta, count -> [ meta, count.text.trim().toInteger() ] }.view()


    emit:
    filter_summary = ch_filter_summary
    versions = ch_versions
}

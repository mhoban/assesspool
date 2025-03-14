include { VCFLIB_VCFFILTER            } from '../../modules/nf-core/vcflib/vcffilter/main'
include { VCFTOOLS as VCFTOOLS_FILTER } from '../../modules/nf-core/vcftools/main'
include { VCFTOOLS as VCFTOOLS_THIN   } from '../../modules/nf-core/vcftools/main'

workflow FILTER {

    take:
    ch_vcf // channel: [ val(meta), [ vcf ] ]

    main:

    ch_versions = Channel.empty()


    def vcff  = params.min_mapping_quality || params.min_mapping_quality_ref || params.min_mapping_ratio ||
                params.max_mapping_ratio || params.read_balance_left || params.read_balance_right ||
                params.quality_depth_ratio || params.mispaired_reads || params.min_pools ||
                params.min_depth || params.max_allele_length || params.min_quality ||
                params.variant_type || params.min_alternate_observations

    if (vcff) {
        VCFLIB_VCFFILTER( ch_vcf.combine(Channel.fromPath('fake.tbi')) )
        ch_versions = ch_versions.mix(VCFLIB_VCFFILTER.out.versions.first())
        VCFLIB_VCFFILTER.out.vcf.set{ ch_vcf2 }
    } else {
        ch_vcf.set{ ch_vcf2 }
    }

    def vcft =  params.max_missing || params.min_minor_allele_count ||params.min_mean_depth ||
            params.max_mean_depth || params.hwe_cutoff

    if (vcft) {
        VCFTOOLS_FILTER ( ch_vcf2, [], [] )
        ch_versions = ch_versions.mix(VCFTOOLS_FILTER.out.versions.first())
        VCFTOOLS_FILTER.out.vcf.set{ ch_vcf_final }
    } else {
        ch_vcf2.set{ ch_vcf_final }
    }

    if (params.thin_snps) {
        VCFTOOLS_THIN(ch_vcf_final)
        ch_versions = ch_versions.mix(VCFTOOLS_THIN.out.versions.first())
        VCFTOOLS_THIN.out.vcf.set{ ch_vcf_final }
    }


    // TODO: visualize filter results? (rmd section 'filter_visualization')


    emit:
    vcf      = ch_vcf_final
    versions = ch_versions                     // channel: [ versions.yml ]
}


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


    // get VCF sample names from header
    VCF_SAMPLES(ch_vcf.combine(Channel.fromPath('fake.tbi')),[],[],[]).output
        .splitText()
        .map { meta, pop -> [ meta.id, pop.trim() ] }
        .groupTuple()
        .set { ch_samples }
    ch_versions = ch_versions.mix(VCF_SAMPLES.out.versions.first())


    // smash VCF sample names into populations
    // if they weren't specified in the samplesheet
    ch_vcf
        .map{ meta, vcf -> [ meta.id, meta, vcf ] }
        .join(ch_samples)
        .map { id, meta, vcf, pops ->
            def pp = meta.populations ?: pops
            def ps = meta.pool_sizes.size() > 1 ? meta.pool_sizes : (1..pops.size()).collect{ meta.pool_sizes[0] }
            [ [ id: meta.id, populations: [pp,ps].transpose().collectEntries(), pop_map: [pops,pp].transpose().collectEntries() ], vcf ]
        }
        .set{ ch_vcf }


    ch_vcf
        .map{ meta, vcf ->
            meta.pop_map
                .collect{ k, v -> "${k} ${v}\n"}
        }
        .flatten()
        .collectFile(name: "sample_remap.txt")
        .ifEmpty(false)
        .set{ ch_remap }

    ch_vcf
        .combine(ch_remap)
        .map{ meta, vcf, samples -> [ meta, vcf, [], samples ]}
        .set{ ch_rename }

    // TODO: use ext.when to control whether we rename or not
    /* e.g.:
    process {
        withName: 'BAR' {
            ext.when = { meta.rename }
        }
    }
    */
    VCF_RENAME(ch_rename,ch_vcf.map{ meta, vcf -> [meta,[]]})
        .vcf
        .set { ch_vcf }
    ch_versions = ch_versions.mix(VCF_RENAME.out.versions.first())

    emit:
    vcf = ch_vcf
    ref = ch_ref

    versions = ch_versions                     // channel: [ versions.yml ]
}


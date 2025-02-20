process GRENEDALF_SYNC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::grenedalf=0.6.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grenedalf:0.6.2--hbefcdb2_1':
        'biocontainers/grenedalf:0.6.2--hbefcdb2_1' }"

    input:
    tuple val(meta), path(vcf), path(index)
    tuple val(meta), path(fasta)
    tuple val(meta), path(fai)

    output:
    tuple val(meta), path("*.sync"), emit: sync
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_arg = fasta ? "--reference-genome-fasta ${fasta}" : ""
    def fai_arg = fai ? "--reference-genome-fai ${fai}" : ""
    """
    grenedalf sync \\
        --threads ${task.cpus} \\
        ${args} \\
        ${fasta_arg} \\
        ${fai_arg} \\
        --file-prefix "${prefix}_" \\
        --vcf-path ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grenedalf: \$(grenedalf --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}_sync.sync

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grenedalf: \$(grenedalf --version)
    END_VERSIONS
    """
}

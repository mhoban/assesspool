process GRENEDALF_SYNC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::grenedalf=0.6.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grenedalf:0.6.2--hbefcdb2_1':
        'biocontainers/grenedalf:0.6.2--hbefcdb2_1' }"

    input:
    tuple val(meta), path(vcf), path(index)
    tuple val(meta), path(sync)
    tuple val(meta), path(fasta)
    tuple val(meta), path(fai)
    tuple val(meta), path(sample_map)

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
    def input_opt = vcf ? '--vcf-path' : (sync ? '--sync-path' : '')
    def input_arg = vcf ?: (sync ?: '')
    def remap_arg = sample_map ? "--rename-samples-list ${sample_map}" : ""
    """
    grenedalf sync \\
        --threads ${task.cpus} \\
        ${args} \\
        ${fasta_arg} \\
        ${fai_arg} \\
        ${remap_arg} \\
        ${input_opt} ${input_arg} \\
        --file-prefix "${prefix}_"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grenedalf: \$(grenedalf --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sync.sync

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grenedalf: \$(grenedalf --version)
    END_VERSIONS
    """
}

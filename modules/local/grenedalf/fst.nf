
process GRENEDALF_FST {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::grenedalf=0.6.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grenedalf:0.6.2--hbefcdb2_1':
        'biocontainers/grenedalf:0.6.2--hbefcdb2_1' }"

    input:
    tuple val(meta), path(sync), path(pool_sizes)
    tuple val(meta), path(fasta)
    tuple val(meta), path(fai)

    output:
    tuple val(meta), path("*fst.csv"), emit: fst
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_arg = fasta ? "--reference-genome-fasta ${fasta}" : ""
    def fai_arg = fai ? "--reference-genome-fai ${fai}" : ""
    """
    grenedalf fst \\
        --file-prefix "${prefix}_" \\
        --threads ${task.cpus} \\
        ${fasta_arg} \\
        ${fai_arg} \\
        ${args} \\
        --sync-path ${sync} \\
        --pool-sizes ${pool_sizes}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grenedalf: \$(grenedalf --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_fst.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grenedalf: \$(grenedalf --version)
    END_VERSIONS
    """
}

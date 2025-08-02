process EXTRACT_SEQUENCES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/mulled-v2-61c59287265e27f2c4589cfc90013ef6c2c6acf1:fb3e48060a8c0e5108b1b60a2ad7e090cfb9eee5-0':
    //     'biocontainers/mulled-v2-61c59287265e27f2c4589cfc90013ef6c2c6acf1:fb3e48060a8c0e5108b1b60a2ad7e090cfb9eee5-0' }"
    container 'fishbotherer/extractr:1.0.0'

    input:
    tuple val(meta), path(fasta), path(fai)
    tuple val(meta), path(fst)
    val(fst_cutoff)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extractsequences.R \\
        $args \\
        --fst-cutoff ${fst_cutoff} \\
        --output ${prefix}.fasta \\
        --index ${fai} \\
        $fasta \\
        $fst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extractsequences: \$(extractsequences --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extractsequences: \$(extractsequences --version)
    END_VERSIONS
    """
}

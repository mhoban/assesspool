process EXTRACT_SEQUENCES {
    tag "$meta.id"
    label 'process_single_mem'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-24c2519e6f79268102473627bccca112fa9f03b6:b565a85a493491c260acda69b780ef6ffa1cda87-0':
        'biocontainers/mulled-v2-24c2519e6f79268102473627bccca112fa9f03b6:b565a85a493491c260acda69b780ef6ffa1cda87-0' }"

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
        R: \$(Rscript -e "cat(paste(R.version[c('major','minor')],collapse='.'))")
        optparse: \$(Rscript -e "cat(paste(packageVersion('optparse')),sep='.')")
        Rsamtools: \$(Rscript -e "cat(paste(packageVersion('Rsamtools')),sep='.')")
        data.table: \$(Rscript -e "cat(paste(packageVersion('data.table')),sep='.')")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(paste(R.version[c('major','minor')],collapse='.'))")
        optparse: \$(Rscript -e "cat(paste(packageVersion('optparse')),sep='.')")
        Rsamtools: \$(Rscript -e "cat(paste(packageVersion('Rsamtools')),sep='.')")
        data.table: \$(Rscript -e "cat(paste(packageVersion('data.table')),sep='.')")
    END_VERSIONS
    """
}

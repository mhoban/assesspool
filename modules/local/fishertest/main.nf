process FISHERTEST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-91e4157032c00e9d9ac47f9837a67abee8f5afc2:7945c7dcdead60dfcbf06b19f4758634d6ad230a-0':
        'biocontainers/mulled-v2-91e4157032c00e9d9ac47f9837a67abee8f5afc2:7945c7dcdead60dfcbf06b19f4758634d6ad230a-0' }"

    input:
    tuple val(meta), val(pools), path(frequency)

    output:
    tuple val(meta), path("*_window_*_fisher.tsv"), emit: fisher
    tuple val(meta), path("*_all_snps_fisher.tsv"), emit: all_snps, optional: true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pool_arg = "--pools '${pools.join(',')}'"
    """
    fisher.R \\
        ${args} \\
        ${pool_arg} \\
        --prefix "${prefix}" \\
        ${frequency}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(paste(R.version[c('major','minor')],collapse='.'))")
        optparse: \$(Rscript -e "cat(paste(packageVersion('optparse')),sep='.')")
        data.table: \$(Rscript -e "cat(paste(packageVersion('data.table')),sep='.')")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_window_stub_fisher.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(paste(R.version[c('major','minor')],collapse='.'))")
        optparse: \$(Rscript -e "cat(paste(packageVersion('optparse')),sep='.')")
        data.table: \$(Rscript -e "cat(paste(packageVersion('data.table')),sep='.')")
    END_VERSIONS
    """
}

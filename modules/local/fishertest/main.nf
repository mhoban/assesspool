process FISHERTEST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-dcb5179d6045c4df41d4319c99d446681271222e:4856497379063722ac7efef973f14ae775bd88ca-0':
        'biocontainers/mulled-v2-dcb5179d6045c4df41d4319c99d446681271222e:4856497379063722ac7efef973f14ae775bd88ca-0' }"

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
        optparse: \$(Rscript -e "cat(paste(packageVersion('dplyr')),sep='.')")
        tibble: \$(Rscript -e "cat(paste(packageVersion('tidyr')),sep='.')")
        janitor: \$(Rscript -e "cat(paste(packageVersion('readr')),sep='.')")
        readr: \$(Rscript -e "cat(paste(packageVersion('janitor')),sep='.')")
        purrr: \$(Rscript -e "cat(paste(packageVersion('purrr')),sep='.')")
        dplyr: \$(Rscript -e "cat(paste(packageVersion('optparse')),sep='.')")
        dplyr: \$(Rscript -e "cat(paste(packageVersion('matrixstats')),sep='.')")
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
        optparse: \$(Rscript -e "cat(paste(packageVersion('dplyr')),sep='.')")
        tibble: \$(Rscript -e "cat(paste(packageVersion('tidyr')),sep='.')")
        janitor: \$(Rscript -e "cat(paste(packageVersion('readr')),sep='.')")
        readr: \$(Rscript -e "cat(paste(packageVersion('janitor')),sep='.')")
        purrr: \$(Rscript -e "cat(paste(packageVersion('purrr')),sep='.')")
        dplyr: \$(Rscript -e "cat(paste(packageVersion('optparse')),sep='.')")
        dplyr: \$(Rscript -e "cat(paste(packageVersion('matrixstats')),sep='.')")
    END_VERSIONS
    """
}

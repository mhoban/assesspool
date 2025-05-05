process JOINFREQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1021c2bc41756fa99bc402f461dad0d1c35358c1:b0c847e4fb89c343b04036e33b2daa19c4152cf5-0':
        'biocontainers/mulled-v2-1021c2bc41756fa99bc402f461dad0d1c35358c1:b0c847e4fb89c343b04036e33b2daa19c4152cf5-0' }"

    input:
    tuple val(meta), path(frequency), path(fst), val(method)

    output:
    tuple val(meta), path("*_fst_frequency.tsv"), emit: fst_freq
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = (task.ext.prefix ?: "${meta.id}") + (method ? "_${method}" : "_unspecified")
    def m = method ? "--method ${method}" : ""
    """
    joinfreq.R \\
        ${args} \\
        ${m} \\
        --prefix "${prefix}" \\
        "${frequency}" \\
        "${fst}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(paste(R.version[c('major','minor')],collapse='.'))")
        dplyr: \$(Rscript -e "cat(paste(packageVersion('dplyr')),sep='.')")
        tidyr: \$(Rscript -e "cat(paste(packageVersion('tidyr')),sep='.')")
        readr: \$(Rscript -e "cat(paste(packageVersion('readr')),sep='.')")
        purrr: \$(Rscript -e "cat(paste(packageVersion('purrr')),sep='.')")
        stringr: \$(Rscript -e "cat(paste(packageVersion('stringr')),sep='.')")
        optparse: \$(Rscript -e "cat(paste(packageVersion('optparse')),sep='.')")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_fst_frequency.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(paste(R.version[c('major','minor')],collapse='.'))")
        dplyr: \$(Rscript -e "cat(paste(packageVersion('dplyr')),sep='.')")
        tidyr: \$(Rscript -e "cat(paste(packageVersion('tidyr')),sep='.')")
        readr: \$(Rscript -e "cat(paste(packageVersion('readr')),sep='.')")
        purrr: \$(Rscript -e "cat(paste(packageVersion('purrr')),sep='.')")
        stringr: \$(Rscript -e "cat(paste(packageVersion('stringr')),sep='.')")
        optparse: \$(Rscript -e "cat(paste(packageVersion('optparse')),sep='.')")
    END_VERSIONS
    """
}

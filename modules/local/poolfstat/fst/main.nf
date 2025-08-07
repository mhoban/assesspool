process POOLFSTAT_FST {
    tag "$meta.id"
    label 'process_medium_low'

    conda "${moduleDir}/environment.yml"
    container 'fishbotherer/poolfstat:1.0.2'
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/mulled-v2-5d425a6e26ab032a4dcba5bc997f43ab9299830a:a1958d454debf3df6f397d04d5f094df68682ef9-0':
    //     'biocontainers/mulled-v2-5d425a6e26ab032a4dcba5bc997f43ab9299830a:a1958d454debf3df6f397d04d5f094df68682ef9-0' }"

    input:
    tuple val(meta), val(pool_map), path(sync)

    output:
    tuple val(meta), path('*.fst')       , emit: fst
    tuple val(meta), path('*global*.tsv'), emit: global
    path('*sliding*.tsv')                , emit: sliding, optional: true
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ps = pool_map.collect{ it.value }.join(',')
    def pn = pool_map.collect{ it.key }.join(',')
    """
    poolfstat.R \\
        ${args} \\
        --pool-sizes ${ps} \\
        --pool-names ${pn} \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${sync}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(paste(R.version[c('major','minor')],collapse='.'))")
        poolfstat: \$(Rscript -e "cat(paste(packageVersion('poolfstat')),sep='.')")
        optparse: \$(Rscript -e "cat(paste(packageVersion('optparse')),sep='.')")
        data.table: \$(Rscript -e "cat(paste(packageVersion('data.table')),sep='.')")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(paste(R.version[c('major','minor')],collapse='.'))")
        poolfstat: \$(Rscript -e "cat(paste(packageVersion('poolfstat')),sep='.')")
        optparse: \$(Rscript -e "cat(paste(packageVersion('optparse')),sep='.')")
        data.table: \$(Rscript -e "cat(paste(packageVersion('data.table')),sep='.')")
    END_VERSIONS
    """
}

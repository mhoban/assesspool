process POOLFSTAT_FST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'quay.io/fishbotherer/poolfstat:1.0.0'

    input:
    tuple val(meta), path(sync), val(pool_map)

    output:
    tuple val(meta), path('*.fst'), path('*global*.tsv'), emit: fst
    path('*sliding*.tsv')                               , emit: sliding, optional: true
    path 'versions.yml'                                , emit: versions

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
    END_VERSIONS
    """
}

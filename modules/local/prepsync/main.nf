process PREPSYNC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'quay.io/fishbotherer/tidy-variantannotation:1.0.0' }"
    container 'quay.io/fishbotherer/tidy-variantannotation:1.0.1'

    input:
    tuple val(meta), path(vcf)
    path(vcf_ref)
    output:
    tuple val(meta), path("combined/*.sync") , optional: true, emit: combined_sync
    tuple val(meta), path("split/*.sync")    , optional: true, emit: split_sync

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prep_sync.R \\
        --prefix ${prefix} \\
        ${args} \\
        ${vcf} \\
        ${vcf_ref}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(paste(R.version[c('major','minor')],collapse='.'))")
    \$(Rscript -e "as.data.frame(installed.packages())[,c(1,ncol(installed.packages()))]" | egrep '^(VariantAnnotation|tidyverse|janitor)' | awk '{print "    " \$1 ": " \$NF}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sync

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(paste(R.version[c('major','minor')],collapse='.'))")
    \$(Rscript -e "as.data.frame(installed.packages())[,c(1,ncol(installed.packages()))]" | egrep '^(VariantAnnotation|tidyverse|janitor)' | awk '{print "    " \$1 ": " \$NF}')
    END_VERSIONS
    """
}

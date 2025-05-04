include { combn } from "./combn"

process POPOOLATION2_FISHERTEST {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::popoolation2=1.201"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popoolation2:1.201--pl5321hdfd78af_0':
        'biocontainers/popoolation2:1.201--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(sync), val(pool_map)

    output:
    tuple val(meta), path("*.fisher"), emit: fisher
    tuple val(meta), path("*.params"), emit: fisher_params
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pools = combn(pool_map.collect{ it.key }, 2).collect { it.sort().join(':') + ".fisher" }.join("\t")
    """
    fisher-test.pl \\
        --input "${sync}" \\
        --output "${sync.BaseName}.fisher" \\
        ${args}

        for fisher in *.fisher; do
            sed -i \$'1i chrom\tmiddle\tsnps\tcovered_fraction\tavg_min_coverage\t${pools}' \$fisher
        done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popoolation2: 1.201
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${sync.BaseName}.fisher"
    touch "${sync.BaseName}.fisher.params"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popoolation2: 1.201
    END_VERSIONS
    """
}

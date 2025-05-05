include { combn } from "./combn"

process POPOOLATION2_FST {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::popoolation2=1.201"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popoolation2:1.201--pl5321hdfd78af_0':
        'biocontainers/popoolation2:1.201--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), val(pool_map), path(sync)

    output:
    tuple val(meta), path("*.fst")   , emit: fst
    tuple val(meta), path("*.params"), emit: fst_params
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ps = pool_map.collect{ it.value }.join(':')
    def pools = combn(pool_map.keySet(), 2).collect { it.sort().join(':') + ".fst" }.join("\t")
    """
    fst-sliding.pl \\
        --input ${sync} \\
        --output "${sync.BaseName}.fst" \\
        --pool-size ${ps} \\
        ${args}

    for fst in *.fst; do
        sed -i \$'1i chrom\tpos\tswindow_size\tcovered_fraction\tavg_min_coverage\t${pools}' \$fst
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
    touch "${sync.BaseName}.fst"
    touch "${sync.BaseName}.fst.params"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popoolation2: 1.201
    END_VERSIONS
    """
}

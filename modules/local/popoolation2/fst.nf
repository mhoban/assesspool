include { combn } from "./combn"

process POPOOLATION2_FST {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::popoolation2=1.201"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popoolation2:1.201--pl5321hdfd78af_0':
        'biocontainers/popoolation2:1.201--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(sync)

    output:
    tuple val(meta), path("*.fst*"), emit: fst
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pool_sizes = meta.populations.collect{ it.value }.join(':')
    def populations = combn(meta.populations.collect{ it.key }, 2).collect { it.join(':') }.join("\t")
    def hdr = task.ext.args2 ?: false
    """
    fst-sliding.pl \\
        --input ${sync} \\
        --output "${sync.BaseName}.fst" \\
        --pool-size ${pool_sizes} \\
        ${args}

    if ${hdr}; then
        for fst in *.fst; do
            sed -i \$'1i chrom\tpos\twindow_size\tcovered_fraction\tavg_min_coverage\t${populations}' \$fst
        done
    fi

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popoolation2: 1.201
    END_VERSIONS
    """
}

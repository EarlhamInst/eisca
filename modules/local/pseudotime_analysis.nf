process PSEUDOTIME_ANALYSIS {
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy_palantir:latest' :
        'docker.io/myeihub/scanpy_palantir:202609' }"

    input:
    path h5ad_annotated
    // path model_file

    output:
    path "pseudotime"
    // path "pseudotime/*.h5ad",  emit: h5ad
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''

    """
    pseudotime_analysis.py \\
        --h5ad ${h5ad_annotated} \\
        --outdir pseudotime \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """



    // stub:
    // """
    // touch combined_matrix.h5ad
    // """
}

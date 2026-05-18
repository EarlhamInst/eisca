process ANNOTATE_CELLS_SCVI {
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scvi_tools_scanpy:latest' :
        'docker.io/myeihub/scvi_tools_scanpy:1.3.3' }"

    input:
    path h5ad_filtered
    // path h5ad_ref

    output:
    path "annotation_scvi"
    path "annotation_scvi/*.h5ad",  emit: h5ad
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    annotate_cells_scvi.py \\
        --h5ad ${h5ad_filtered} \\
        --outdir annotation_scvi \\
        $args \\


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

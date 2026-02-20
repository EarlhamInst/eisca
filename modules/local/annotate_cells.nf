process ANNOTATE_CELLS {
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/celltypist:1.6.3--pyhdfd78af_0' :
        'docker.io/myeihub/celltypist_ps:1.6.3' }"

    input:
    path h5ad_filtered
    // path model_file

    output:
    path "annotation/celltypist"
    path "annotation/celltypist/*.h5ad",  emit: h5ad
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    
    // Conditionally build the model argument
    def model_arg = ''
    if (params.ctmodel) {
        // Ensure you handle the case where params.ctmodel might be a path object or string
        model_arg = "--model_file ${params.ctmodel}"
    }

    """
    annotate_cells.py \\
        --h5ad ${h5ad_filtered} \\
        --outdir annotation/celltypist \\
        ${model_arg} \\
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

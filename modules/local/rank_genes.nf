process RANK_GENES {
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy-scripts:1.9.301--pyhdfd78af_0' :
        'biocontainers/scanpy-scripts:1.9.301--pyhdfd78af_0' }"

    input:
    path(h5ad_filtered)
    // path samplesheet

    output:
    path "dea/${subfolder}"
    // path "clustering/*.h5ad",  emit: h5ad
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (args.contains('--celltype_col')) {
        subfolder = 'compare_ct'
    } else if (args.contains('--groupby group')) {
        subfolder = 'compare'
    } else if (args.contains('--combine')) {
        subfolder = 'markers_cb'
    } else {
        subfolder = 'markers'
    }

    """
    rank_genes.py \\
        --h5ad ${h5ad_filtered} \\
        --outdir dea/${subfolder} \\
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

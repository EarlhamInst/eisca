process CELLCELL_COMMUNICATION {
    label 'process_high'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/davele23/cellchat_sc_rst:latest' :
        'docker.io/myeihub/cellchat_sc_rst:2.1.2' }"

    input:
    path(countmtx)
    path(metadata)
    path(gene_ids)
    path(cell_ids)
    // path samplesheet

    output:
    path "cellchat"
    // path "clustering/*.h5ad",  emit: h5ad
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cellcell_communication.R \\
        --count ${countmtx} \\
        --metadata ${metadata} \\
        --gids ${gene_ids} \\
        --cids ${cell_ids} \\
        --outdir cellchat \\
        $args \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellchat: 2.1.2
    END_VERSIONS        
    """



    // stub:
    // """
    // touch combined_matrix.h5ad
    // """
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_eisca_pipeline'


include { KALLISTO_BUSTOOLS                  } from '../subworkflows/local/kallisto_bustools'
include { SCRNASEQ_ALEVIN                    } from '../subworkflows/local/alevin'
include { STARSOLO                           } from '../subworkflows/local/starsolo'
include { MTX_CONVERSION                     } from "../subworkflows/local/mtx_conversion"
include { GTF_GENE_FILTER                    } from '../modules/local/gtf_gene_filter'
//include { EMPTYDROPS_CELL_CALLING            } from '../modules/local/emptydrops'
include { GUNZIP as GUNZIP_FASTA             } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF               } from '../modules/nf-core/gunzip/main'
include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'
include { getGenomeAttribute                 } from '../subworkflows/local/utils_nfcore_eisca_pipeline'

include { QC_CELL_FILTER                    } from '../modules/local/qc_cell_filter'
include { CLUSTERING_ANALYSIS               } from '../modules/local/clustering_analysis'
include { ANNOTATE_CELLS                    } from '../modules/local/annotate_cells'
include { ANNOTATE_CELLS_SCVI               } from '../modules/local/annotate_cells_scvi'
include { DEA_SCVI                          } from '../modules/local/dea_scvi'
include { TRAIN_CT_MODEL                    } from '../modules/local/train_ct_model'
include { RANK_GENES                        } from '../modules/local/rank_genes'
include { TRIMGALORE                        } from '../modules/nf-core/trimgalore/main'
include { H5AD_TO_MTX_META                  } from '../modules/local/h5ad_to_mtx_meta'
include { CELLCELL_COMMUNICATION            } from '../modules/local/cellcell_communication'
// include { MAKE_REPORT                       } from '../modules/local/make_report'
// include { GET_PARAMS                       } from '../modules/local/get_params'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// analyses = params.analyses.split(',').toList()
// skip_analyses = params.skip? params.skip.split(',').toList() : []

workflow EISCA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    // samplesheet - this is passed to the MTX conversion functions to add metadata to the AnnData objects.
    ch_input = file(params.input)

    if (params.run_analyses.contains('primary')){
        protocol_config = Utils.getProtocol(workflow, log, params.aligner, params.protocol)
        if (protocol_config['protocol'] == 'auto' && params.aligner !in ["cellranger", "cellrangerarc", "cellrangermulti"]) {
            error "Only cellranger supports `protocol = 'auto'`. Please specify the protocol manually!"
        }

        ch_genome_fasta = params.fasta ? file(params.fasta, checkIfExists: true) : ( params.genome ? file( getGenomeAttribute('fasta'), checkIfExists: true ) : [] )
        ch_gtf          = params.gtf   ? file(params.gtf  , checkIfExists: true) : ( params.genome ? file( getGenomeAttribute('gtf')  , checkIfExists: true ) : [] )

        // general input and params
        ch_transcript_fasta = params.transcript_fasta ? file(params.transcript_fasta): []
        //ch_motifs = params.motifs ? file(params.motifs) : []
        //ch_cellrangerarc_config = params.cellrangerarc_config ? file(params.cellrangerarc_config) : []
        ch_txp2gene = params.txp2gene ? file(params.txp2gene) : []
        ch_multiqc_files = Channel.empty()
        if (params.barcode_whitelist) {
            ch_barcode_whitelist = file(params.barcode_whitelist)
        } else if (protocol_config.containsKey("whitelist")) {
            ch_barcode_whitelist = file("$projectDir/${protocol_config['whitelist']}")
        } else {
            ch_barcode_whitelist = []
        }
    

        //kallisto params
        ch_kallisto_index = params.kallisto_index ? file(params.kallisto_index) : []
        kb_workflow = params.kb_workflow
        kb_t1c = params.kb_t1c ? file(params.kb_t1c) : []
        kb_t2c = params.kb_t2c ? file(params.kb_t2c) : []

        //salmon params
        ch_salmon_index = params.salmon_index ? file(params.salmon_index) : []

        //star params
        ch_star_index = params.star_index ? file(params.star_index) : []
        star_feature = params.star_feature
    }


    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_mtx_matrices = Channel.empty()

    
    //===================================== Primary anaysis stage =====================================

    if (params.run_analyses.contains('primary')){
    
        // MODULE: Run FastQC
        if (!params.skip_analyses.contains('fastqc')) {
            FASTQC (
                ch_samplesheet
            )
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        }


        // Uncompress genome fasta file if required
        if (params.fasta) {
            if (params.fasta.endsWith('.gz')) {
                ch_genome_fasta    = GUNZIP_FASTA ( [ [:], file(params.fasta) ] ).gunzip.map { it[1] }
                ch_versions        = ch_versions.mix(GUNZIP_FASTA.out.versions)
            } else {
                ch_genome_fasta = Channel.value( file(params.fasta) )
            }
        }

        //
        // Uncompress GTF annotation file or create from GFF3 if required
        //
        if (params.gtf) {
            if (params.gtf.endsWith('.gz')) {
                ch_gtf      = GUNZIP_GTF ( [ [:], file(params.gtf) ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
            } else {
                ch_gtf = Channel.value( file(params.gtf) )
            }
        }

        // filter gtf
        ch_filter_gtf = ch_gtf ? GTF_GENE_FILTER ( ch_genome_fasta, ch_gtf ).gtf : []


        // Trimming reads
        if (params.trim) {
            TRIMGALORE (ch_samplesheet)
            ch_samplesheet = TRIMGALORE.out.reads
            ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
        }


        // Run kallisto bustools pipeline
        if (params.aligner == "kallisto") {
            KALLISTO_BUSTOOLS(
                ch_genome_fasta,
                ch_filter_gtf,
                ch_kallisto_index,
                ch_txp2gene,
                kb_t1c,
                kb_t2c,
                protocol_config['protocol'],
                kb_workflow,
                ch_samplesheet
            )
            ch_versions = ch_versions.mix(KALLISTO_BUSTOOLS.out.ch_versions)
            ch_mtx_matrices = ch_mtx_matrices.mix(KALLISTO_BUSTOOLS.out.raw_counts, KALLISTO_BUSTOOLS.out.filtered_counts)
            if (!params.protocol.contains('smartseq')){  //avoid spaced path generated from KALLISTO_BUSTOOLS
                ch_txp2gene = KALLISTO_BUSTOOLS.out.txp2gene
            }
        }

        // Run salmon alevin pipeline
        if (params.aligner == "alevin") {
            SCRNASEQ_ALEVIN(
                ch_genome_fasta,
                ch_filter_gtf,
                ch_transcript_fasta,
                ch_salmon_index,
                ch_txp2gene,
                ch_barcode_whitelist,
                protocol_config['protocol'],
                ch_samplesheet
            )
            ch_versions = ch_versions.mix(SCRNASEQ_ALEVIN.out.ch_versions)
            ch_multiqc_files = ch_multiqc_files.mix(SCRNASEQ_ALEVIN.out.alevin_results.map{ meta, it -> it })
            ch_mtx_matrices = ch_mtx_matrices.mix(SCRNASEQ_ALEVIN.out.alevin_results)
        }

        // Run STARSolo pipeline
        if (params.aligner == "star") {
            STARSOLO(
                ch_genome_fasta,
                ch_filter_gtf,
                ch_star_index,
                protocol_config['protocol'],
                ch_barcode_whitelist,
                ch_samplesheet,
                star_feature,
                protocol_config.get('extra_args', ""),
            )
            ch_versions = ch_versions.mix(STARSOLO.out.ch_versions)
            ch_mtx_matrices = ch_mtx_matrices.mix(STARSOLO.out.raw_counts, STARSOLO.out.filtered_counts)
            ch_star_index = STARSOLO.out.star_index
            ch_multiqc_files = ch_multiqc_files.mix(STARSOLO.out.for_multiqc)
        }


        // Run mtx to h5ad conversion subworkflow
        MTX_CONVERSION (
            ch_mtx_matrices,
            ch_input,
            ch_txp2gene,
            ch_star_index
        )

        //Add Versions from MTX Conversion workflow too
        ch_versions.mix(MTX_CONVERSION.out.ch_versions)

    }



    //===================================== Secondary anaysis stage =====================================

    if (params.run_analyses.contains('secondary')){
    
        // MODULE: Run QC and cell filtering
        ch_h5ad = Channel.empty()
        if(params.run_analyses.contains('primary')){
            ch_h5ad = MTX_CONVERSION.out.h5ad
        }else if(params.h5ad){
            ch_h5ad = Channel.fromPath(params.h5ad)
        }else if(params.aligner){
            path = [
                'kallisto': "${params.outdir}/kallisto/mtx_conversions/combined_*_matrix.h5ad",
                'alevin': "${params.outdir}/alevin/mtx_conversions/combined_*_matrix.h5ad",
                'star': "${params.outdir}/star/mtx_conversions/combined_raw_matrix.h5ad"
            ].get(params.aligner)
            path2 = "${params.outdir}/qc_cell_filter/adata_filtered_normalized.h5ad"
            if(params.skip_analyses.contains('qccellfilter') && (new File(path2).exists())){
                ch_h5ad = Channel.fromPath(path2)
            }else{
                ch_h5ad = Channel.fromPath(path)
            } 
        }else{
            log.warn("For this analysis, please specify an h5ad file either by setting --aligner for the " +
            "h5ad file generated by the aligner or by setting --h5ad for an existing h5ad file.")
            return
        }

        if (!params.skip_analyses.contains('qccellfilter')) {
            QC_CELL_FILTER (
                ch_h5ad,
                Channel.fromPath(params.input)
                // MTX_CONVERSION.out.h5ad
            )
            // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(QC_CELL_FILTER.out.versions)
            ch_h5ad = QC_CELL_FILTER.out.h5ad         
        }
        
        if (!params.skip_analyses.contains('clustering')) {
            CLUSTERING_ANALYSIS (
                ch_h5ad
            )
            ch_versions = ch_versions.mix(CLUSTERING_ANALYSIS.out.versions)
        }   
        
    }



    //===================================== Tertiary anaysis stage =====================================

    if (params.run_analyses.any{ it in [
        'tertiary', 
        'annotation', 
        'annotation_scvi', 
        'dea', 
        'dea_scvi', 
        'cellchat'] }){
    
        // Get input h5ad file
        ch_h5ad = Channel.empty()
        if(params.run_analyses.contains('secondary')){
            if (!params.skip_analyses.contains('clustering')) {
                ch_h5ad = CLUSTERING_ANALYSIS.out.h5ad
            }else {
                ch_h5ad = QC_CELL_FILTER.out.h5ad
            }
        }else if(params.h5ad){
            ch_h5ad = Channel.fromPath(params.h5ad)
        }else{
            path1 = "${params.outdir}/clustering/adata_clustering.h5ad"
            path2 = "${params.outdir}/qc_cell_filter/adata_filtered_normalized.h5ad"
            path3 = "${params.outdir}/annotation/adata_annotation.h5ad"
            path4 = "${params.outdir}/annotation_scvi/adata_annotation.h5ad"
            if(params.run_analyses.any{ it in ['dea', 'dea_scvi', 'cellchat'] }){
                if(new File(path3).exists()){
                    ch_h5ad = Channel.fromPath(path3)
                }else if(new File(path4).exists()){
                    ch_h5ad = Channel.fromPath(path4)
                }     
            }else if(new File(path1).exists()){
                ch_h5ad = Channel.fromPath(path1)
            }else if(new File(path2).exists()){
                ch_h5ad = Channel.fromPath(path2)
            }
        }
        ch_h5ad.ifEmpty {
            log.warn("For this analysis, h5ad file can be found in secondary analysis, please specify an h5ad file by setting --h5ad.")
            return            
        }

        if (params.run_analyses.any{it=='tertiary' || it=='annotation'} && !params.skip_analyses.contains('annotation')) {
            // ch_ctmodel = params.ctmodel? Channel.fromPath(params.ctmodel) : []
            ANNOTATE_CELLS (
                ch_h5ad,
                // ch_ctmodel
                // MTX_CONVERSION.out.h5ad
            )
            // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(ANNOTATE_CELLS.out.versions)
            ch_h5ad = ANNOTATE_CELLS.out.h5ad      
        }

        if (params.run_analyses.any{it=='tertiary' || it=='annotation_scvi'} && 
            !params.skip_analyses.contains('annotation_scvi') && params.args_annotation_scvi) {
            // ch_h5ad_ref = params.h5ad_ref? Channel.fromPath(params.h5ad_ref) : []
            ANNOTATE_CELLS_SCVI (
                ch_h5ad,
                // ch_h5ad_ref,
                // ch_ctmodel
                // MTX_CONVERSION.out.h5ad
            )
            // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(ANNOTATE_CELLS_SCVI.out.versions)
            ch_h5ad = ANNOTATE_CELLS_SCVI.out.h5ad      
        }
        
        if (params.run_analyses.any{it=='tertiary' || it=='dea'} && !params.skip_analyses.contains('dea')) {
            RANK_GENES (
                ch_h5ad,
            )
            ch_versions = ch_versions.mix(RANK_GENES.out.versions)
        }

        if (params.run_analyses.any{it=='tertiary' || it=='dea_scvi'} && 
            !params.skip_analyses.contains('dea_scvi') && params.args_dea_scvi) {
            // ch_h5ad_ref = params.h5ad_ref? Channel.fromPath(params.h5ad_ref) : []
            DEA_SCVI (
                ch_h5ad,
                // ch_h5ad_ref,
                // ch_ctmodel
                // MTX_CONVERSION.out.h5ad
            )
            // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(DEA_SCVI.out.versions)
        }

        if (params.run_analyses.any{it=='tertiary' || it=='cellchat'} && !params.skip_analyses.contains('cellchat')) {
            H5AD_TO_MTX_META (
                ch_h5ad,
            )
            ch_versions = ch_versions.mix(H5AD_TO_MTX_META.out.versions)
            CELLCELL_COMMUNICATION (
                H5AD_TO_MTX_META.out.h5ad_mtx,
                H5AD_TO_MTX_META.out.h5ad_meta,
                H5AD_TO_MTX_META.out.h5ad_genes,
                H5AD_TO_MTX_META.out.h5ad_cells,
            )
            ch_versions = ch_versions.mix(CELLCELL_COMMUNICATION.out.versions)
        }
        
    }

    // train cell-type models for CellTypist
    if (params.run_analyses.contains('ctmodel') && !params.skip_analyses.contains('ctmodel')) {
            ch_h5ad = Channel.fromPath(params.h5ad)
            TRAIN_CT_MODEL (
                ch_h5ad,
                // MTX_CONVERSION.out.h5ad
            )
            // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(TRAIN_CT_MODEL.out.versions)
            // ch_h5ad = ANNOTATE_CELLS.out.h5ad      
    }





    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    // done = true
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

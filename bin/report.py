#!/usr/bin/env python

import dominate.tags as html
import ezcharts as ezc
from ezcharts.components.reports.labs import LabsReport, LabsNavigation, ILabsNavigationClasses
from ezcharts.layout.snippets import DataTable, Grid, Tabs
from ezcharts.components.theme import LAB_head_resources
from ezcharts.components.ezchart import EZChart

import argparse
from pathlib import Path
from report_util import *
import util
import json
import pandas as pd
import sys


logger = util.get_named_logger('REPORT')

report_title = 'EI Single Cell Analysis Report'
workflow_name = 'eisca'


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Quality control before and after cell filtering",
        # epilog="python count_reads_from_bam.py --bam file.bam --bed file.bed --json output.json",
    )
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--results",
        metavar="RESULTS_DIR",
        type=Path,
        help="Results directory.",
        required=True,
    )
    # parser.add_argument(
    #     "--stats", nargs='+',
    #     help="Fastcat per-read stats, ordered as per entries in --metadata.")
    # parser.add_argument(
    #     "--images", nargs='+',
    #     help="Sample directories containing various images to put in report")
    # parser.add_argument(
    #     "--survival",
    #     help="Read survival data in TSV format")
    parser.add_argument(
        "--params",
        metavar="FILE_PARAMS",
        type=Path,
        help="Workflow params json file",
        required=True,
    )
    parser.add_argument(
        "--versions",
        metavar="FILE_VERSIONS",
        type=Path,
        help="Workflow versions file",
        required=True,
    )
    parser.add_argument(
        "--logo",
        metavar="FILE_LOGO",
        type=Path,
        help="Logo image file",
        required=True,
    )    
    # parser.add_argument(
    #     "--umap_dirs", nargs='+',
    #     help="Sample directories containing umap and gene expression files")
    # parser.add_argument(
    #     "--umap_genes", help="File containing list of genes to annnotate UMAPs")
    # parser.add_argument(
    #     "--metadata", default='metadata.json', required=True,
    #     help="sample metadata")
    parser.add_argument(
        "--wf-version", default='unknown',
        help="version of the executed workflow")
    parser.add_argument(
        "--samplesheet",
        metavar="SAMPLESHEET",
        type=Path,
        help="Input samplesheet file.",
        required=True,
    )           
    return parser.parse_args(argv)


def main(argv=None):
    logger.info('Building report')
    args = parse_args(argv)


    if not args.params.is_file():
        logger.error(f"The given input file {args.params} was not found!")
        sys.exit(2)

    if not args.versions.is_file():
        logger.error(f"The given input file {args.versions} was not found!")
        sys.exit(2)

    report = LabsReport(
        report_title, workflow_name,
        args.params, args.versions, args.wf_version,
        head_resources=[*LAB_head_resources],
    )
    
    report.nav.getElementsByTagName('a')[0].clear()
    navdiv = report.nav.getElementsByTagName('div')[0]
    link = html.a(href='https://www.earlham.ac.uk/', cls=ILabsNavigationClasses().logo)
    with link:
        # with open(Path('images/EI_logo.png'), 'rb') as f:
        with open(Path(args.logo), 'rb') as f:
            b64img = base64.b64encode(f.read()).decode()
            html.img(src=f'data:image/png;base64,{b64img}', width=120)
    navdiv[0] = link

    report.banner.clear()
    report.footer.clear()
    report.intro_content.add(EIBanner(report_title, workflow_name))
    report.footer.add(EILabsAddendum(workflow_name, args.wf_version))

    samplesheet = pd.read_csv(args.samplesheet)
    # batch = 'group' if 'group' in samplesheet.columns else 'sample'
    sample = 'plate' if 'plate' in samplesheet.columns else  'sample'
    batch = 'sample'
    if 'group' in samplesheet.columns:
        batch = 'group'
    elif 'plate' in samplesheet.columns:
        batch = 'plate'  
    # Nbatch = len(samplesheet[batch].unique())
    # samples = samplesheet['sample'].unique()

    path_quant_qc = Path(args.results, 'qc_cell_filter')
    path_quant_qc_raw = Path(path_quant_qc, 'raw_counts')
    # path_quant_qc_scatter = Path(path_quant_qc, 'scatter')
    # path_quant_qc_violin = Path(path_quant_qc, 'violin')
    path_cell_filtering = Path(path_quant_qc, 'cell_filtering')
    # path_cell_filtering_dist = Path(path_cell_filtering, 'distribution')
    path_clustering = Path(args.results, 'clustering')
    path_annotation = Path(args.results, 'annotation')
    path_annotation_scvi = Path(args.results, 'annotation_scvi')
    path_dea = Path(args.results, 'dea')
    path_dea_scvi = Path(args.results, 'dea_scvi')
    path_cellchat = Path(args.results, 'cellchat')

    # print(path_quant_qc) #tst

    if path_quant_qc.exists():
        summary = pd.read_csv(Path(path_quant_qc, 'sample_summary.csv')).set_index(f"{sample.capitalize()} ID")
        Nsample = summary.shape[0]
        with report.add_section('Single cell summary', 'Summary'):
            html.p("""This section gives an overall summary of the single-cell count matrix for 
                   each sample. The statistics include the total number of cells with at least 
                   one gene expressed, the total number of genes expressed in at least one cell, 
                   the median number of genes per cell, and the median percentage of counts in 
                   mitochondrial genes (pct-mt).""")
            DataTable.from_pandas(summary)   

        with report.add_section('Quantification QC', 'QC'):
            html.p("""This section presents the QC plots of the raw count matrix generated during 
                   the quantification step. These plots provide insight into the quality of the 
                   experiments and guide the filtering of low-quality cells.""")
            html.p("""The following scatter plot shows the relationship between total 
                   read counts and the number of genes, with the percentage of counts in 
                   mitochondrial genes indicated by color.""")
            plots_from_image_files(path_quant_qc_raw, meta='sample', widths=['800'], suffix=['scatter*.png'])
            html.p("""The following knee plot shows the relationship between cell barcode rank and total 
                   number of UMIs per barcode. The inflection (knee) point where the UMI frequency rapidly 
                   drops can be good cutoff to remove empty droplets.""")
            plots_from_image_files(path_quant_qc_raw, meta='sample', widths=['500'], suffix=['knee_plot.png'])
            html.p("""The following violin plots display the distribution of cells based on the number of 
                   genes, total counts, and the percentage of counts in mitochondrial genes.""")
            plots_from_image_files(path_quant_qc_raw, meta='sample', ncol=3, suffix=['violin*.png'])
    else:
        logger.info('Skipping Quantification QC')

    if path_cell_filtering.exists():
        summary = pd.read_csv(Path(path_cell_filtering, 'sample_summary_filtered.csv')).set_index(f"{sample.capitalize()} ID")
        with report.add_section('Cell filtering', 'Cell filtering'):
            html.p("""This section presents the statistics and QC plots after cell filtering process. 
                   The filtering process includes hard thresholds for minimum genes, minimum counts, 
                   minimum cells and the percentage of counts in mitochondrial genes. Additionally, 
                   users can set quantile limits on the number of genes. Then doublet detection is 
                   performed using Scrublet.""")
            html.p("""The following table shows summary statistics, with percentages in brackets 
                   indicating the comparison to the raw counts.""")
            DataTable.from_pandas(summary)
            html.p("""The following violin plots display the distribution of cells based on the number of 
                   genes, total counts, and the percentage of counts in mitochondrial genes after filtering.""")
            plots_from_image_files(path_cell_filtering, meta='sample', ncol=3, suffix=['violin*.png'])            
            # html.p("""The following plots show the distribution KDE curves before and after filtering 
            #        for the number of genes, total counts, and the percentage of counts in mitochondrial genes.""")            
            # plots_from_image_files(path_cell_filtering, meta='sample', ncol=2, suffix=['dist*.png'])
            html.p("""The following plots show the UMAP plots
                   for the number of genes, total counts, and the percentage of counts in mitochondrial genes.""")                        
            plots_from_image_files(path_cell_filtering, meta='sample', suffix=['umap_total*.png'])
            if util.check_file(f"{path_cell_filtering}/sample_*", 'umap_doublet.png'):
                html.p("""The following plots show the UMAP plots for the predicted doublets and doublet scores.""")                        
                plots_from_image_files(path_cell_filtering, meta='sample', suffix=['umap_doublet.png'])
            show_analysis_parameters(f"{path_quant_qc}/parameters.json")          
    else:
        logger.info('Skipping Cell filtering')


    if path_clustering.exists():
        if util.check_file(f"{path_clustering}/sample_*", ''):
            batch = 'sample'
        elif util.check_file(f"{path_clustering}/group_*", ''):
            batch = 'group'
        Nbatch = len(samplesheet[batch].unique())      
        with report.add_section('Clustering analysis', 'Clustering'):
            html.p(f"""This section shows clustering UMAP plots for each {batch}. The clustering 
                   was performed using Leiden graph-clustering method. The resolution parameter 
                   was set for different values to get different number of clusters which 
                   could match to biologically-meaningful cell types.""")
            plots_from_image_files(path_clustering, meta=batch, ncol=2)
            html.p(f"""The following plot shows a stacked bar chart that presents the proportions of clusters 
                   across {batch}s, calculated for each resolution value. The plot illustrates the distribution 
                   profiles of predicted clusters between {batch}s.""")                        
            plots_from_image_files(path_clustering, meta='resolution', widths=[str(min(Nbatch*280, 1200))])
            show_analysis_parameters(f"{path_clustering}/parameters.json")                   
    else:
        logger.info('Skipping clustering analysis')        


    if path_annotation.exists() or path_annotation_scvi.exists():
        with report.add_section('Cell-type annotation', 'Annotation'):
            if path_annotation.exists():
                if util.check_file(f"{path_annotation}/sample_*", ''):
                    batch = 'sample'
                elif util.check_file(f"{path_annotation}/group_*", ''):
                    batch = 'group'
                Nbatch = len(samplesheet[batch].unique())               
                html.p(f"""This section presents cell-type annotation results using CellTypist which is an 
                automated tool for cell type annotation based on pre-trained models, capable of accurately 
                classifying different cell types and subtypes.""")
                html.p("""The following UMAP plots show the predicted cell-type clusters and the mapped 
                confidence scores of the cells.""")            
                plots_from_image_files(path_annotation, suffix=['umap_cell_type.png'], meta=batch, widths=['1200'])
                plots_from_image_files(path_annotation, suffix=['umap_conf_score.png'], meta=batch, widths=['600'])
                html.p(f"""The following plot shows a stacked bar chart that presents the proportions 
                    of cell-type clusters across {batch}s. The plot illustrates the distribution 
                    profiles of predicted cell-type clusters between {batch}s.""")                   
                plots_from_image_files(path_annotation, suffix=['prop_*.png'], widths=[str(min(Nbatch*500, 1200))])
                show_analysis_parameters(f"{path_annotation}/parameters.json")
                if path_annotation_scvi.exists(): html.hr(style="border: 1px solid grey;")

            if path_annotation_scvi.exists():
                if util.check_file(f"{path_annotation_scvi}/sample_*", ''):
                    batch = 'sample'
                elif util.check_file(f"{path_annotation_scvi}/group_*", ''):
                    batch = 'group'
                Nbatch = len(samplesheet[batch].unique())               
                html.p(f"""This section presents cell-type annotation results using scvi-tools, which predicts the 
                       cell types of the query dataset based on a scANVI model trained on the reference latent 
                       space and transferred to the query data.""")
                html.p("""The following UMAP plots show the predicted cell-type clusters and the mapped confidence scores 
                       which are probabilities associated with that predicted cell types.""")            
                plots_from_image_files(path_annotation_scvi, suffix=['umap_cell_type.png'], meta=batch, widths=['1200'])
                plots_from_image_files(path_annotation_scvi, suffix=['umap_conf_score.png'], meta=batch, widths=['600'])
                html.p(f"""The following plot shows a stacked bar chart that presents the proportions 
                    of cell-type clusters across {batch}s. The plot illustrates the distribution 
                    profiles of predicted cell-type clusters between {batch}s.""")                   
                plots_from_image_files(path_annotation_scvi, suffix=['prop_*.png'], widths=[str(min(Nbatch*500, 1200))])
                show_analysis_parameters(f"{path_annotation_scvi}/parameters.json")                            
    else:
        logger.info('Skipping cell-type annotation')   


    if path_dea.exists() or path_dea_scvi.exists():
        with report.add_section('Differential expression analysis', 'DEA'):
            if path_dea.exists():
                if util.check_file(f"{path_dea}/sample_*", ''):
                    batch = 'sample'
                elif util.check_file(f"{path_dea}/group_*", ''):
                    batch = 'group'      
            
                html.p("""This section presents the results of the differentially expression analysis using Scanpy's 
                    rank_genes_groups function. These results allow users to identify marker genes by comparing 
                    the ranked genes of one cluster against all others, as well as to explore differentially 
                    expressed genes between two conditions.""")

                # showing plots for DEA between conditions for all cells
                if util.check_file(f"{path_dea}", '*.png'):
                    html.p("""The following plots show differentially expressed genes between the two conditions.""")                        
                    plots_from_image_files(path_dea, suffix=['plot_genes_*.png'])
                    plots_from_image_files(path_dea, suffix=['dotplot_genes_*.png'])

                # showing plots for DEA between conditions for each celltype
                if util.check_file(f"{path_dea}/celltype_*", '*.png'):
                    html.p("""The following plots show differentially expressed genes between two conditions across various cell types/clusters.""")                        
                    plots_from_image_files(path_dea, meta='celltype', suffix=['plot_genes_*.png'])
                    plots_from_image_files(path_dea, meta='celltype', suffix=['dotplot_genes_*.png'])

                # showing plots for one cluster vs rest for each sample/group
                if util.check_file(f"{path_dea}/{batch}_*", '*.png'):
                    html.p(f"""The following plots display the ranking of genes for one of the cell clusters against the rest of the clusters across {batch}s.""")                        
                    plots_from_image_files(path_dea, meta=batch, suffix=['plot_genes_*.png'])
                    plots_from_image_files(path_dea, meta=batch, suffix=['dotplot_genes_*.png'])

                show_analysis_parameters(f"{path_dea}/parameters.json")
                if path_annotation_scvi.exists(): html.hr(style="border: 1px solid grey;")

            if path_dea_scvi.exists():
                if util.check_file(f"{path_dea_scvi}/sample_*", ''):
                    batch = 'sample'
                elif util.check_file(f"{path_dea_scvi}/group_*", ''):
                    batch = 'group'      
            
                html.p("""This section presents the results of the differentially expression analysis using scvi-tools's 
                       differential_expression function by estimating the posterior distribution of the log fold-change (LFC) 
                       between subpopulations. These results allow users to identify marker genes by comparing the ranked 
                       genes of one cluster against all others, as well as to explore differentially expressed genes between 
                       two conditions.""")

                # showing plots for DEA between conditions for all cells
                if util.check_file(f"{path_dea_scvi}", '*.png'):
                    html.p("""The following plots show differentially expressed genes between the two conditions.""")                        
                    plots_from_image_files(path_dea_scvi, suffix=['dotplot_*.png'])
                    plots_from_image_files(path_dea_scvi, suffix=['heatmap_*.png'])

                # showing plots for DEA between conditions for each celltype
                if util.check_file(f"{path_dea_scvi}/celltype_*", '*.png'):
                    html.p("""The following plots show differentially expressed genes between two conditions across selected subpopulations.""")                        
                    plots_from_image_files(path_dea_scvi, meta='celltype', suffix=['dotplot_*.png'])
                    plots_from_image_files(path_dea_scvi, meta='celltype', suffix=['heatmap_*.png'])

                # showing plots for one cluster vs rest for each sample/group
                if util.check_file(f"{path_dea_scvi}/{batch}_*", '*.png'):
                    html.p(f"""The following plots show differentially expressed genes between subpopulations across {batch}s.""")                        
                    plots_from_image_files(path_dea_scvi, meta=batch, suffix=['dotplot_*.png'])
                    plots_from_image_files(path_dea_scvi, meta=batch, suffix=['heatmap_*.png'])

                show_analysis_parameters(f"{path_dea_scvi}/parameters.json")
    else:
        logger.info('Skipping differential expression analysis')


    if path_cellchat.exists():
        if util.check_file(f"{path_cellchat}/sample_*", ''):
            batch = 'sample'
        elif util.check_file(f"{path_cellchat}/group_*", ''):
            batch = 'group'        
        with report.add_section('Cell-cell communication analysis', 'CellChat'):
            html.p("""This section presents the results of the cell–cell communication analysis using CellChat. These 
                   results enable users to identify, analyze, and visualize intercellular communication networks between 
                   cell groups from scRNA-seq data, based on the CellChat database of ligand–receptor interactions.""")

            # showing the aggregated cell-cell communication network
            if util.check_file(f"{path_cellchat}/{batch}_*", 'aggregated_network_all.png'):
                html.p(f"""The following plots show the aggregated cell-cell communication network, either by counting the 
                       number of links or by summarizing communication probabilities (weights). These circle plots provide 
                       an overview of cell-cell interactions: the left plot shows the number of interactions, while the right 
                       plot displays the total interaction strength between cell groups. Thicker edge lines indicate stronger 
                       signals, and the circle sizes are proportional to the number of cells in each group.""")                        
                plots_from_image_files(path_cellchat, meta=batch, ncol=2, suffix=['aggregated_network_all*.png'])

            # showing the aggregated cell-cell communication network for each signaling pathway
            if util.check_file(f"{path_cellchat}/{batch}_*", 'pathway_network_circle_*.png'):
                html.p(f"""The following plots show the aggregated (weighted) cell-cell communication networks for 
                       individual signalling pathways with significant ligand-receptor interactions.""")                        
                plots_from_image_files(path_cellchat, meta=batch, ncol=3, suffix=['pathway_network_circle_*.png'])

            # showing the aggregated cell-cell communication network for each signaling pathway (heatmap)
            # if util.check_file(f"{path_cellchat}/{batch}_*", 'pathway_network_heatmap_*.png'):
            #     html.p(f"""The following heatmaps show interaction strength in the cell-cell communication network between 
            #            pairs of cell groups, with the y-axis representing sender cells and the x-axis representing receiver 
            #            cells. The top bar plot shows the column-wise sums of the heatmap values, while the right bar plot 
            #            shows the row-wise sums.""")                        
            #     plots_from_image_files(path_cellchat, meta=batch, ncol=2, suffix=['pathway_network_heatmap_*.png'])

            # showing significant interactions (L-R pairs) from each cell group to other cell groups
            if util.check_file(f"{path_cellchat}/{batch}_*", 'cellcell_LR_bubble_*.png'):
                html.p(f"""The following bubble plots show all significant interactions (ligand-receptor pairs) from each 
                       cell group to others. The color and size of each dot represent the communication probability and 
                       p-value, respectively. This visualization shows how cell-cell communication is mediated by multiple 
                       ligand-receptor pairs.""")                        
                plots_from_image_files(path_cellchat, meta=batch, ncol=2, suffix=['cellcell_LR_bubble_*.png'])

            show_analysis_parameters(f"{path_cellchat}/parameters.json")                 
    else:
        logger.info('Skipping cell-cell communication analysis')


    report.write(args.report)
    logger.info('Report writing finished')
    

if __name__ == "__main__":
    sys.exit(main())

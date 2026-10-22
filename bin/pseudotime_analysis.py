#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'
import scanpy as sc
import pandas as pd
import numpy as np
import palantir
from matplotlib import pyplot as plt
import argparse
import sys, json
from pathlib import Path
import util

logger = util.get_named_logger('PSEUDOTIME_ANALYSIS')


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    
    parser = argparse.ArgumentParser(
        description="Perform pseudotime and trajectory analysis",
    )
    parser.add_argument(
        "--h5ad",
        metavar="FILE_H5AD",
        type=Path,
        help="Input anndata data file.",
        required=True,
    )
    parser.add_argument(
        "--outdir",
        metavar="OUT_DIR",
        type=Path,
        help="Output directory.",
        required=True,
    )
    parser.add_argument(
        "--n_comp",
        type=int,
        help="Set the number of diffusion components.",
        default=10,
    )
    parser.add_argument(
        "--knn",
        type=int,
        help="Set the number of nearest neighbors for diffusion map.",
        default=30,
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0,
        help="Normalization parameter for the diffusion operator.",
    )
    parser.add_argument(
        "--celltype_col",
        default="cell_type",
        help="Specify a column used to define cell-types or clusters.",
    )
    parser.add_argument(
        "--start_celltype",
        default=None,
        help="Specify a cell-type or cluster for selecting a start cell.",
    )
    parser.add_argument(
        "--end_celltypes",
        default=None,
        help="Specify a list of cell-types for selecting terminal states, e.g. 'celltype1,celltype2'.",
    )
    parser.add_argument(
        "--n_waypoints",
        type=int,
        help="Set the number of waypoints for pseudotime analysis.",
        default=1200,
    )
    parser.add_argument(
        "--q",
        type=float,
        default=0.01,
        help="Quantile threshold for fate assignment for selecting branch cells.",
    )
    parser.add_argument(
        "--eps",
        type=float,
        default=0.01,
        help="Tolerance offset from Quantile fate assignment for selecting branch cells.",
    )
    parser.add_argument(
        "--smoothness",
        type=float,
        default=1.0,
        help="Set the smoothness of the trajectories plot of multiple branches.",
    )
    parser.add_argument(
        "--genes",
        default=None,
        help="Specify a list of genes for gene trend analysis, e.g. 'gene1,gene2'.",
    )
    parser.add_argument(
        "--dm1max",
        help="Select the start cell with the maximum value of diffusion component 1.",
        action='store_true',
    )
    parser.add_argument(
        "--meta",
        default='auto',
        choices=['auto', 'sample', 'group', 'plate'],
        help="Choose a metadata column as the batch for clustering",
    )
    parser.add_argument(
        "--fontsize",
        type=int,
        help="Set font size for plots.",
        default=12,
    )
    parser.add_argument(
        "--pdf",
        help="Whether to generate figure files in PDF format.",
        action='store_true',
    )           
    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)

    if not args.h5ad.is_file():
        logger.error(f"The given input file {args.h5ad} was not found!")
        sys.exit(2)

    plt.rcParams.update({
        "font.size": args.fontsize,
        # "figure.dpi": 100,
        # "axes.titlesize": 'medium',
        # "axes.labelsize": 'small',
        # "xtick.labelsize": 'small',
        # "ytick.labelsize": 'small',
        # "legend.fontsize": 'small',
    })

    util.check_and_create_folder(args.outdir)
    path_outdir = Path(args.outdir)
    util.check_and_create_folder(path_outdir)

    adata = sc.read_h5ad(args.h5ad)
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()

    if args.meta == 'auto':
        # batch = 'group' if hasattr(adata.obs, 'group') else 'sample'
        batch = 'sample'
        if 'group' in adata.obs.columns:
            batch = 'group'
        elif 'plate' in adata.obs.columns:
            batch = 'plate'  
    else:
        batch = args.meta

    if args.celltype_col not in adata.obs.columns:
        logger.error(f"The given cell type column {args.celltype_col} was not found in the data!")
        sys.exit(2)

    # calculate diffusion maps
    dm_res = palantir.utils.run_diffusion_maps(adata, n_components=args.n_comp, knn=args.knn, alpha=args.alpha)
    ms_data = palantir.utils.determine_multiscale_space(adata)
    # uses MAGIC to impute the data
    imputed_X = palantir.utils.run_magic_imputation(adata)

    # find start cell
    if args.start_celltype:
        if args.start_celltype not in adata.obs[args.celltype_col].unique():
            logger.error(f"The given start cell-type {args.start_celltype} was not found in the column {args.celltype_col}!")
            sys.exit(2)
        start_cell = palantir.utils.early_cell(adata, args.start_celltype, args.celltype_col)
    else:
        dm1 = ms_data.iloc[:, 0]
        if args.dm1max:
            start_cell = dm1.idxmax()
        else:
            start_cell = dm1.idxmin()

    # find terminal states
    if args.end_celltypes:
        end_celltypes = args.end_celltypes.split(",")
        for celltype in end_celltypes:
            if celltype not in adata.obs[args.celltype_col].unique():
                logger.error(f"The given end cell-type {celltype} was not found in the column {args.celltype_col}!")
                sys.exit(2)
        terminal_states = palantir.utils.find_terminal_states(adata, end_celltypes, celltype_column=args.celltype_col)
    else:
        terminal_states = palantir.utils.find_terminal_states(adata, adata.obs[args.celltype_col], celltype_column=args.celltype_col)

    # run palantir to compute pseudotime and branch probabilities
    pr_res = palantir.core.run_palantir(adata, start_cell, num_waypoints=args.n_waypoints, terminal_states=terminal_states)
    # select cells associated with a specific branch of the pseudotime trajectory
    masks = palantir.presults.select_branch_cells(adata, q=args.q, eps=args.eps)


    # plot pseudotime, entropy and branch probabilities 
    with plt.rc_context():
        palantir.plot.plot_palantir_results(adata, s=3)
        plt.savefig(Path(path_outdir, f"palantir_results_all.png"), bbox_inches="tight")    
        if args.pdf:
            plt.savefig(Path(path_outdir, f"palantir_results_all.pdf"), bbox_inches="tight") 

    for sid in sorted(adata.obs[batch].unique()):
        adata_s = adata[adata.obs[batch]==sid]   
        path_outdir_s = Path(path_outdir, f"{batch}_{sid}")
        util.check_and_create_folder(path_outdir_s)
        with plt.rc_context({"figure.dpi": 100}):
            palantir.plot.plot_palantir_results(adata_s, s=3)

            plt.savefig(Path(path_outdir_s, f"palantir_results_{sid}.png"), bbox_inches="tight")    
            if args.pdf:
                plt.savefig(Path(path_outdir_s, f"palantir_results_{sid}.pdf"), bbox_inches="tight")  

    # plot trajectories
    with plt.rc_context():
        palantir.plot.plot_trajectories(adata, pseudotime_interval=(0, .9), smoothness=args.smoothness,
            legend_kwargs={
                "bbox_to_anchor": (1.3, 1),
                "loc": "upper left",          
                "frameon": False              
            }
        )
        plt.savefig(Path(path_outdir, f"trajectory_combined.png"), bbox_inches="tight")    
        if args.pdf:
            plt.savefig(Path(path_outdir, f"trajectory_combined.pdf"), bbox_inches="tight")         

    for branch in terminal_states:
        with plt.rc_context():
            palantir.plot.plot_trajectory(adata, branch)
            plt.savefig(Path(path_outdir, f"trajectory_branch_{branch}.png"), bbox_inches="tight")    
            if args.pdf:
                plt.savefig(Path(path_outdir, f"trajectory_branch_{branch}.pdf"), bbox_inches="tight") 

    for sid in sorted(adata.obs[batch].unique()):
        adata_s = adata[adata.obs[batch]==sid]   
        path_outdir_s = Path(path_outdir, f"{batch}_{sid}")
        util.check_and_create_folder(path_outdir_s)
        with plt.rc_context():
            palantir.plot.plot_trajectories(adata_s, pseudotime_interval=(0, .9), smoothness=args.smoothness,
                legend_kwargs={
                    "bbox_to_anchor": (1.3, 1),
                    "loc": "upper left",          
                    "frameon": False              
                }
            )
            plt.savefig(Path(path_outdir_s, f"trajectory_combined_{sid}.png"), bbox_inches="tight")    
            if args.pdf:
                plt.savefig(Path(path_outdir_s, f"trajectory_combined_{sid}.pdf"), bbox_inches="tight")         

        for branch in terminal_states:
            with plt.rc_context():
                palantir.plot.plot_trajectory(adata_s, branch)
                plt.savefig(Path(path_outdir_s, f"trajectory_branch_{branch}_{sid}.png"), bbox_inches="tight")    
                if args.pdf:
                    plt.savefig(Path(path_outdir_s, f"trajectory_branch_{branch}_{sid}.pdf"), bbox_inches="tight") 


    # compute gene trends along trajectory branches using the MAGIC-imputed data
    gene_trends = palantir.presults.compute_gene_trends(adata, expression_key="MAGIC_imputed_data")
    genes = []
    if args.genes:
        genes = args.genes.split(",")
    else:
        for branch in terminal_states:
            genes += (
                gene_trends[branch]["trends"]
                .std(axis=1)
                .sort_values(ascending=False)
                .head(2)
                .index.tolist())

    # plot gene trends along trajectory branches
    for i in range(len(genes)):
        with plt.rc_context():
            palantir.plot.plot_gene_trends(adata, genes=genes[i:i+1])
            plt.savefig(Path(path_outdir, f"gene_trends_{genes[i]}.png"), bbox_inches="tight")    
            if args.pdf:
                plt.savefig(Path(path_outdir, f"gene_trends_{genes[i]}.pdf"), bbox_inches="tight")            

    for sid in sorted(adata.obs[batch].unique()):
        adata_s = adata[adata.obs[batch]==sid]   
        path_outdir_s = Path(path_outdir, f"{batch}_{sid}")
        util.check_and_create_folder(path_outdir_s)

        try:
            gene_trends_s = palantir.presults.compute_gene_trends(adata_s, expression_key="MAGIC_imputed_data")
        except Exception as e:
            logger.error(f"Error occurred while computing gene trends for sample {sid}: {e}")
            continue

        for i in range(len(genes)):
            with plt.rc_context():
                palantir.plot.plot_gene_trends(adata_s, genes=genes[i:i+1])
                plt.savefig(Path(path_outdir_s, f"gene_trends_{sid}_{genes[i]}.png"), bbox_inches="tight")    
                if args.pdf:
                    plt.savefig(Path(path_outdir_s, f"gene_trends_{sid}_{genes[i]}.pdf"), bbox_inches="tight")  


    # save the AnnData into a h5ad file
    adata.write_h5ad(Path(path_outdir, 'adata_palantir.h5ad'), compression="gzip")


    # save analysis parameters into a json file
    with open(Path(path_outdir, 'parameters.json'), 'w') as file:
        params = {}
        params.update({"--h5ad": str(args.h5ad)})         
        params.update({"--n_comp": args.n_comp})        
        params.update({"--knn": args.knn})        
        params.update({"--alpha": args.alpha})        
        params.update({"--celltype_col": args.celltype_col})               
        if args.start_celltype: params.update({"--start_celltype": args.start_celltype})
        if args.end_celltypes: params.update({"--end_celltypes": args.end_celltypes})
        params.update({"--n_waypoints": args.n_waypoints})   
        params.update({"--q": args.q})   
        params.update({"--eps": args.eps})   
        if args.genes: params.update({"--genes": args.genes})
        if args.dm1max: params.update({"--dm1max": args.dm1max})
        params.update({"--meta": args.meta})
        json.dump(params, file, indent=4)



if __name__ == "__main__":
    sys.exit(main())

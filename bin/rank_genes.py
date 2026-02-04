#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import io
import anndata
from matplotlib import pyplot as plt
import argparse
import sys, json
from pathlib import Path
import util

logger = util.get_named_logger('DEA')


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Perform cell clustering and plot UMAPs of clustering",
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
        "--groupby",
        default='leiden',
        help="The key of the observations grouping to consider.",
    )
    parser.add_argument(
        "--groups",
        default='all',
        help="Specify a subset of groups, e.g. 'group1,group2'.",
    )
    parser.add_argument(
        "--reference",
        default='rest',
        help="If spcecify a group name, compare with respect to this group.",
    )    
    parser.add_argument(
        "--method",
        default='t-test',
        choices=['t-test', 'wilcoxon', 'logreg', 't-test_overestim_var'],
        help="Choose a test method for differential expression anlaysis.",
    )
    parser.add_argument(
        "--n_genes",
        type=int,
        default=20,
        help="Number of genes to show in plots",
    )
    parser.add_argument(
        "--meta",
        default='auto',
        choices=['auto', 'sample', 'group', 'plate'],
        help="Choose a metadata column as the batch for clustering.",
    )
    parser.add_argument(
        "--celltype_col",
        default=None,
        help="Spcecify a column used to define cell-types for DEA between groups.",
    )
    parser.add_argument(
        "--celltypes",
        default=None,
        help="Spcecify a list cell-types for DEA between groups, e.g. 'celltype1,celltype2'.",
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
        # "axes.titlesize": 'medium',
        # "axes.labelsize": 'small',
        # "xtick.labelsize": 'small',
        # "ytick.labelsize": 'small',
        # "legend.fontsize": 'small',
    })

    util.check_and_create_folder(args.outdir)
    path_analysis = Path(args.outdir)
    util.check_and_create_folder(path_analysis)

    adata = anndata.read_h5ad(args.h5ad)

    if not adata.uns.get('log1p'): # to fix issue in scanpy function
        adata.uns['log1p'] = {'base': None}

    groupby = args.groupby
    if not hasattr(adata.obs, groupby):
        cols = [col for col in adata.obs.columns if col.startswith(groupby)]
        if cols:
            groupby == cols[0]
        else:
            logger.error(f"Please specify a observation column for grouping!")
            sys.exit(2)

    groups = args.groups.split(',') if args.groups!='all' else None

    if args.meta == 'auto':
        # batch = 'group' if hasattr(adata.obs, 'group') else 'sample'
        batch = 'sample'
        if hasattr(adata.obs, 'group'):
            batch = 'group'
        elif hasattr(adata.obs, 'plate'):
            batch = 'plate'  
    else:
        batch = args.meta

    
    # differential expression analysis
    if groupby == 'group': # between conditions
        if groups == None:
            groups = list(adata.obs['group'].unique())
            groups.remove(args.reference)

        if args.celltype_col: # DEA between conditions for each celltype
            celltypes = args.celltypes.split(',') if args.celltypes else sorted(adata.obs[args.celltype_col].unique()) 
            for celltype in celltypes:
                adata_s = adata[adata.obs[args.celltype_col]==celltype]   
                path_analysis_s = Path(path_analysis, f"celltype_{celltype}".replace(' ', '_').replace('/', '_'))
                util.check_and_create_folder(path_analysis_s)

                sc.tl.rank_genes_groups(
                    adata_s, 
                    groupby, 
                    method=args.method, 
                    groups=groups, 
                    reference=args.reference,
                )
                with plt.rc_context():
                    sc.pl.rank_genes_groups(
                        adata_s, 
                        n_genes=args.n_genes, 
                        sharey=True,
                        groups=groups,
                        fontsize=13,
                    )
                    plt.savefig(Path(path_analysis_s, f"plot_genes_group_{args.reference}.png"), bbox_inches="tight")
                    if args.pdf:
                        plt.savefig(Path(path_analysis_s, f"plot_genes_group_{args.reference}.pdf"), bbox_inches="tight")

                with plt.rc_context():
                    sc.pl.rank_genes_groups_dotplot(
                        adata_s, 
                        n_genes=args.n_genes, 
                        groups=groups,
                    )
                    plt.savefig(Path(path_analysis_s, f"dotplot_genes_group_{args.reference}.png"), bbox_inches="tight")
                    if args.pdf:
                        plt.savefig(Path(path_analysis_s, f"dotplot_genes_group_{args.reference}.pdf"), bbox_inches="tight")

                for gid in groups:
                    sc.get.rank_genes_groups_df(adata_s, group=gid).to_csv(
                        Path(path_analysis_s, f'dea_group_{gid}_vs_{args.reference}.csv'), 
                        index=False,
                    )                            
        else: # DEA between conditions for all cells
            sc.tl.rank_genes_groups(
                adata, 
                groupby, 
                method=args.method, 
                groups=groups, 
                reference=args.reference,
            )
            with plt.rc_context():
                sc.pl.rank_genes_groups(
                    adata, 
                    n_genes=args.n_genes, 
                    sharey=True,
                    groups=groups,
                    fontsize=13,
                )
                plt.savefig(Path(path_analysis, f"plot_genes_group_{args.reference}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_analysis, f"plot_genes_group_{args.reference}.pdf"), bbox_inches="tight")

            with plt.rc_context():
                sc.pl.rank_genes_groups_dotplot(
                    adata, 
                    n_genes=args.n_genes, 
                    groups=groups,
                )
                plt.savefig(Path(path_analysis, f"dotplot_genes_group_{args.reference}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_analysis, f"dotplot_genes_group_{args.reference}.pdf"), bbox_inches="tight")

            for gid in groups:
                sc.get.rank_genes_groups_df(adata, group=gid).to_csv(
                    Path(path_analysis, f'dea_group_{gid}_vs_{args.reference}.csv'), 
                    index=False,
                )
    else: # one cluster vs rest for each sample/group
        for sid in sorted(adata.obs[batch].unique()):
            adata_s = adata[adata.obs[batch]==sid]   
            path_analysis_s = Path(path_analysis, f"{batch}_{sid}")
            util.check_and_create_folder(path_analysis_s)

            # filter out groups which only have one cell
            adata_s = adata_s[adata_s.obs[groupby].astype(str).map(adata_s.obs[groupby].value_counts()) > 1]

            sc.tl.rank_genes_groups(
                adata_s, 
                groupby, 
                method=args.method, 
                groups=groups if groups else 'all', 
                reference=args.reference,
            )
            with plt.rc_context():
                sc.pl.rank_genes_groups(
                    adata_s, 
                    n_genes=args.n_genes, 
                    sharey=True,
                    groups=groups,
                    fontsize=13,
                )
                plt.savefig(Path(path_analysis_s, f"plot_genes_{groupby}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_analysis_s, f"plot_genes_{groupby}.pdf"), bbox_inches="tight")

            with plt.rc_context():
                sc.pl.rank_genes_groups_dotplot(
                    adata_s, 
                    n_genes=args.n_genes, 
                    groups=groups,
                )
                plt.savefig(Path(path_analysis_s, f"dotplot_genes_{groupby}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_analysis_s, f"dotplot_genes_{groupby}.pdf"), bbox_inches="tight")

            for gid in groups if groups else sorted(adata_s.obs[groupby].unique()):
                sc.get.rank_genes_groups_df(adata_s, group=gid).to_csv(
                    Path(path_analysis_s, f'dea_{groupby}_{gid}_vs_{args.reference}.csv'), 
                    index=False,
                )


    # save analysis parameters into a json file
    with open(Path(path_analysis, 'parameters.json'), 'w') as file:
        params = {}
        params.update({"--h5ad": str(args.h5ad)})        
        params.update({"--groupby": args.groupby})
        params.update({"--groups": args.groups})
        params.update({"--reference": args.reference})
        params.update({"--method": args.method})
        params.update({"--n_genes": args.n_genes})
        params.update({"--meta": args.meta})
        if args.celltype_col:
            params.update({"--celltype_col": args.celltype_col}) 
            params.update({"--celltypes": args.celltypes}) 
        json.dump(params, file, indent=4)



if __name__ == "__main__":
    sys.exit(main())

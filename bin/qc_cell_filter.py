#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import scrublet as scr
from scipy import stats
import anndata
from matplotlib import pyplot as plt
import argparse
import sys, json
from pathlib import Path
import util
# from .util import get_named_logger,

logger = util.get_named_logger('QC_CELL_FILTER')



def parse_args(argv=None):
    """Define and immediately parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Quality control before and after cell filtering",
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
        "--samplesheet",
        metavar="SAMPLESHEET",
        type=Path,
        help="Path to samplesheet file.",
    )    
    parser.add_argument(
        "--min_genes",
        type=util.intlist,
        help="Filter cells by minimum number of genes.",
        default=[100],
    )
    parser.add_argument(
        "--max_genes",
        type=util.intlist,
        help="Filter cells by maximum number of genes.",
        default=[0],
    )
    parser.add_argument(
        "--min_counts",
        type=util.intlist,
        help="Filter cells by minimum number of counts.",
        default=[1],
    )
    parser.add_argument(
        "--max_counts",
        type=util.intlist,
        help="Filter cells by maximum number of counts.",
        default=[0],
    )       
    parser.add_argument(
        "--min_cells",
        type=util.intlist,
        help="Filter genes by number of cells expressed.",
        default=[3],
    )
    parser.add_argument(
        "--min_gcounts",
        type=util.intlist,
        help="Filter genes by minimum counts expressed.",
        default=[0],
    )
    parser.add_argument(
        "--pct_mt",
        type=int,
        help="Filter genes by the maximum percentage of mitochondrial counts.",
        default=20,
    )
    parser.add_argument(
        "--doublet_rate",
        type=float,
        help="The expected fraction of transcriptomes that are doublets.",
        default=0.1,
    )        
    parser.add_argument(
        "--quantile_upper",
        type=float,
        help="Filter genes by upper limit of quantile on number of genes.",
        default=1, # 0.98
    )
    parser.add_argument(
        "--quantile_lower",
        type=float,
        help="Filter genes by lower limit of quantile on number of genes.",
        default=0, # 0.02
    )
    parser.add_argument(
        "--find_doublets",
        help="Whether to perform doublets prediction.",
        action='store_true',
    )
    parser.add_argument(
        "--iqr_coef",
        type=float,
        help="Remove outliers which larger than iqr_coef*IQR in total_counts.",
        default=2,
    )
    parser.add_argument(
        "--mt",
        default='MT-',
        help="The prefix of mitochondrial gene IDs")
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

    if not args.samplesheet.is_file():
        logger.error(f"The given input file {args.samplesheet} was not found!")
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
    # path_quant_qc = Path(args.outdir, 'quant_qc')
    path_quant_qc = Path(args.outdir)
    path_quant_qc_raw = Path(path_quant_qc, 'raw_counts')
    # path_quant_qc_scatter = Path(path_quant_qc, 'scatter')
    # path_quant_qc_violin = Path(path_quant_qc, 'violin')
    
    path_cell_filtering = Path(args.outdir, 'cell_filtering')
    # path_cell_filtering_dist = Path(path_cell_filtering, 'distribution')
    util.check_and_create_folder(path_quant_qc)
    util.check_and_create_folder(path_quant_qc_raw)
    # util.check_and_create_folder(path_quant_qc_scatter)
    # util.check_and_create_folder(path_quant_qc_violin)
    util.check_and_create_folder(path_cell_filtering)

    samplesheet = pd.read_csv(args.samplesheet)
    if 'group' in samplesheet.columns: # check if all samples assigned a group
        if sum(samplesheet['group'].isna()) > 0:
            logger.error(f"In the 'group' column, all samples must be assigned to a group.")
            sys.exit(1)

    adata_raw = sc.read_h5ad(args.h5ad)
    # adata_raw.obs['sample'] = [x.removesuffix('_raw') for x in adata_raw.obs['sample']]

    summary = []
    summary_filtered = []
    adatas = {}
    sample = 'plate' if hasattr(adata_raw.obs, 'plate') else 'sample'
    # for sid in adata_raw.obs['sample'].unique():
    for sid in samplesheet[sample].unique():
        adatas.update({sid: adata_raw[adata_raw.obs[sample]==sid]})     

    # for sid, adata_s in adatas.items():
    for s, (sid, adata_s) in enumerate(adatas.items()):

        # QC on raw counts
        # mitochondrial genes
        adata_s.var["mt"] = adata_s.var_names.str.startswith(args.mt)
        # ribosomal genes
        adata_s.var["ribo"] = adata_s.var_names.str.startswith(("RPS", "RPL"))
        # hemoglobin genes
        adata_s.var["hb"] = adata_s.var_names.str.contains("^HB[^(P)]")

        sc.pp.calculate_qc_metrics(
            adata_s, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
        )

        # obs_raw = adata_s.obs.copy()

        # create summary csv file for all samples
        n_cells_raw = adata_s.obs[adata_s.obs['n_genes_by_counts']>0].shape[0]
        n_genes_raw = adata_s.var[adata_s.var['n_cells_by_counts']>0].shape[0]
        summary += [{
            f"{sample.capitalize()} ID": sid,
            'Number of cells': n_cells_raw,
            'Number of genes': n_genes_raw,
            'Median genes per cell': np.nanmedian(adata_s.obs['n_genes_by_counts']),
            'Median of pct-mt':  np.nanmedian(adata_s.obs['pct_counts_mt']),
        }]

        # scatter plots on total_counts against n_genes_by_counts
        path_sample = Path(path_quant_qc_raw, f"sample_{sid}")
        util.check_and_create_folder(path_sample)

        with plt.rc_context():
            sc.pl.scatter(adata_s, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
            plt.savefig(Path(path_sample, 'scatter_total_counts_genes.png'), bbox_inches="tight")
            if args.pdf:
                plt.savefig(Path(path_sample, 'scatter_total_counts_genes.pdf'), bbox_inches="tight")
        
        with plt.rc_context():
            counts = np.sort(adata_s.obs["total_counts"])[::-1]
            ranks = np.arange(1, len(counts) + 1)
            plt.figure(figsize=(6, 5))
            plt.plot(ranks, counts, lw=2)
            plt.xscale("log") #base-10 logarithm
            plt.yscale("log")
            plt.xlabel("Cell barcode rank")
            plt.ylabel("Total counts (UMI frequency)")
            plt.title("Knee plot")
            plt.grid(True, which="both", lw=0.2)
            plt.savefig(Path(path_sample, 'knee_plot.png'), bbox_inches="tight")
            if args.pdf:
                plt.savefig(Path(path_sample, 'knee_plot.pdf'), bbox_inches="tight")


        # violin plots for n_genes_by_counts, total_counts, pct_counts_mt
        path_sample = Path(path_quant_qc_raw, f"sample_{sid}")
        util.check_and_create_folder(path_sample)
        for i, qc in enumerate(["n_genes_by_counts", "total_counts", "pct_counts_mt"]):
            with plt.rc_context():
                sc.pl.violin(
                    adata_s,
                    [qc],
                    jitter=0.4,
                    show=False
                )
                plt.savefig(Path(path_sample, f'violin{i}_{qc}.png'), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_sample, f'violin{i}_{qc}.pdf'), bbox_inches="tight")

        # due to bug in old scanpy if using multi_panel=True
        # with plt.rc_context():
        #     sc.pl.violin(
        #         adata_s,
        #         ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        #         jitter=0.4,
        #         multi_panel=True,
        #         show=False
        #     )
        #     plt.savefig(Path(path_sample, 'violin_total_counts_genes_mt.png'), bbox_inches="tight")


        # Cell filtering
        min_genes_s = args.min_genes[s] if s < len(args.min_genes) else args.min_genes[-1]
        min_counts_s = args.min_counts[s] if s < len(args.min_counts) else args.min_counts[-1]
        min_cells_s = args.min_cells[s] if s < len(args.min_cells) else args.min_cells[-1]
        min_gcounts_s = args.min_gcounts[s] if s < len(args.min_gcounts) else args.min_gcounts[-1]
        max_counts_s = args.max_counts[s] if s < len(args.max_counts) else args.max_counts[-1]
        max_genes_s = args.max_genes[s] if s < len(args.max_genes) else args.max_genes[-1]

        sc.pp.filter_cells(adata_s, min_genes=min_genes_s)
        sc.pp.filter_cells(adata_s, min_counts=min_counts_s)
        sc.pp.filter_genes(adata_s, min_cells=min_cells_s)
        if min_gcounts_s > 0:
            sc.pp.filter_genes(adata_s, min_counts=min_gcounts_s)
        if max_counts_s > 0:
            sc.pp.filter_cells(adata_s, max_counts=max_counts_s)
        if max_genes_s > 0:
            sc.pp.filter_cells(adata_s, max_genes=max_genes_s)

        sc.pp.calculate_qc_metrics(adata_s, inplace=True)

        if args.iqr_coef > 0:
            q1 = np.percentile(adata_s.obs.total_counts.values, 25)
            q3 = np.percentile(adata_s.obs.total_counts.values, 75)
            upper_fence = q3 + args.iqr_coef*(q3 - q1)
            sc.pp.filter_cells(adata_s, max_counts=upper_fence)
            sc.pp.calculate_qc_metrics(adata_s, inplace=True)
        # if args.iqr_coef > 0:        
        #     q1 = np.percentile(adata_s.obs.n_genes_by_counts.values, 25)
        #     q3 = np.percentile(adata_s.obs.n_genes_by_counts.values, 75)
        #     upper_fence = q3 + args.iqr_coef*(q3 - q1)
        #     sc.pp.filter_cells(adata_s, max_genes=upper_fence)

        if args.quantile_upper < 1:
            upper_lim = np.quantile(adata_s.obs.n_genes_by_counts.values, args.quantile_upper)
            adata_s = adata_s[adata_s.obs.n_genes_by_counts.values < upper_lim]
        if args.quantile_lower > 0:
            lower_lim = np.quantile(adata_s.obs.n_genes_by_counts.values, args.quantile_lower)
            adata_s = adata_s[adata_s.obs.n_genes_by_counts.values > lower_lim]
        
        adata_s = adata_s[adata_s.obs.pct_counts_mt < args.pct_mt]


        # Doublet detection
        if args.find_doublets:
            try:
                scrub = scr.Scrublet(adata_s.X, expected_doublet_rate=args.doublet_rate)
                doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
                adata_s.obs['doublet_score'] = doublet_scores
                adata_s.obs['predicted_doublet'] = predicted_doublets
                adata_s.obs['predicted_doublet_1'] = [int(x) for x in predicted_doublets] # for ploting due to bugs in scanpy
            except ValueError as e:
                logger.warning(f"Failed to perform doublet detection due to error: {e}")


        # create summary csv file for all samples after filtering
        obs_s = adata_s.obs[adata_s.obs['n_genes_by_counts']>0].copy()
        if 'predicted_doublet' in adata_s.obs.columns:
            obs_s = obs_s[~obs_s['predicted_doublet']]
        n_cells = obs_s.shape[0]
        n_genes = adata_s.var[adata_s.var['n_cells_by_counts']>0].shape[0]
        summary_filtered += [{
            f"{sample.capitalize()} ID": sid,
            'Number of cells': f"{n_cells} ({n_cells/n_cells_raw:.0%})",
            'Number of genes': f"{n_genes} ({n_genes/n_genes_raw:.0%})",
            'Median genes per cell': np.nanmedian(obs_s['n_genes_by_counts']),
            'Median of pct-mt':  np.nanmedian(obs_s['pct_counts_mt']),
        }]


        # distributions of n_genes_by_counts, total_counts and pct_counts_mt before and after filtering
        # obs1 = obs_raw[obs_raw.n_genes_by_counts>0]   
        # obs2 = adata_s.obs.copy() 
        # if 'predicted_doublet' in obs2.columns:
        #     obs2 = obs2[~obs2['predicted_doublet']]
        path_sample = Path(path_cell_filtering, f"sample_{sid}")
        util.check_and_create_folder(path_sample)

        # violin plots for n_genes_by_counts, total_counts, pct_counts_mt
        for i, qc in enumerate(["n_genes_by_counts", "total_counts", "pct_counts_mt"]):
            with plt.rc_context():
                sc.pl.violin(
                    adata_s[~adata_s.obs['predicted_doublet']] if hasattr(adata_s.obs, 'predicted_doublet') else adata_s,
                    [qc],
                    jitter=0.4,
                    show=False
                )
                plt.savefig(Path(path_sample, f'violin{i}_{qc}.png'), bbox_inches="tight")        
                if args.pdf:
                    plt.savefig(Path(path_sample, f'violin{i}_{qc}.pdf'), bbox_inches="tight")        

        # hist_df = pd.concat([obs1[['n_genes_by_counts']], obs2[['n_genes_by_counts']]], axis=1)
        # hist_df.columns = ['Before filtering', 'After filtering']
        # plot = sns.displot(data=hist_df, kind="kde",  fill=True, palette=sns.color_palette('bright'), aspect=1.2)
        # plot.set(xlabel='n_genes_by_counts')
        # plot.ax.set_xlim(right=max(obs2['n_genes_by_counts']))
        # plot.add_legend(loc="upper right")
        # plot.savefig(Path(path_sample, 'dist1_n_genes_by_counts.png'))

        # hist_df = pd.concat([obs1[['total_counts']], obs2[['total_counts']]], axis=1)
        # hist_df.columns = ['Before filtering', 'After filtering']
        # plot = sns.displot(data=hist_df, kind="kde",  fill=True, palette=sns.color_palette('bright'), aspect=1.2)
        # plot.set(xlabel='total_counts')
        # plot.ax.set_xlim(right=max(obs2['total_counts']))
        # plot.add_legend(loc="upper right")
        # plot.savefig(Path(path_sample, 'dist2_total_counts.png'))

        # hist_df = pd.concat([obs1[['pct_counts_mt']], obs2[['pct_counts_mt']]], axis=1)
        # hist_df.columns = ['Before filtering', 'After filtering']
        # plot = sns.displot(data=hist_df, kind="kde",  fill=True, palette=sns.color_palette('bright'), aspect=1.2)
        # plot.set(xlabel='pct_counts_mt')
        # plot.ax.set_xlim(right=max(obs2['pct_counts_mt']))
        # plot.add_legend(loc="upper right")
        # plot.savefig(Path(path_sample, 'dist3_pct_counts_mt.png'))


        # save a filtered and normalized h5ad file for each sample
        # adata_s.write_h5ad(Path(path_quant_qc, f'adata_filtered_normalized_{sid}.h5ad'))

        adatas[sid] = adata_s


    summary = pd.DataFrame(summary, index=adatas.keys())
    summary.to_csv(Path(path_quant_qc, 'sample_summary.csv'), index=False)
    summary_filtered = pd.DataFrame(summary_filtered, index=adatas.keys()) 
    summary_filtered.to_csv(Path(path_cell_filtering, 'sample_summary_filtered.csv'), index=False)

    # combine all samples' anndata into one anndata
    adata = anndata.concat(adatas, label="sample", index_unique="_", join="outer", merge="unique")
    adata.layers["counts"] = adata.X.copy()  # preserve raw counts

    # Normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Feature selection and dimensionality reduction
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    with plt.rc_context():
        sc.pl.highly_variable_genes(adata, show=False)
        plt.savefig(Path(path_cell_filtering, 'highly_variable_genes.png'), bbox_inches="tight")
        if args.pdf:
            plt.savefig(Path(path_cell_filtering, 'highly_variable_genes.pdf'), bbox_inches="tight")
    sc.tl.pca(adata)

    # QC after cell filtering
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    with plt.rc_context():
        sc.pl.umap(
            adata,
            color="sample",
            size=2,
            show=False
        )
        plt.savefig(Path(path_cell_filtering, 'umap_samples.png'), bbox_inches="tight")
        if args.pdf:
            plt.savefig(Path(path_cell_filtering, 'umap_samples.pdf'), bbox_inches="tight")


    for sid in adata.obs[sample].unique():
        adata_s = adata[adata.obs[sample]==sid]
        path_cell_filtering_s = Path(path_cell_filtering, f"sample_{sid}")
        util.check_and_create_folder(path_cell_filtering_s)

        if 'predicted_doublet' in adata_s.obs.columns:
            with plt.rc_context():
                sc.pl.umap(
                    adata_s,
                    color=["predicted_doublet_1", "doublet_score"],
                    wspace=0.3, # increase horizontal space between panels
                    size=3,
                    show=False
                )
                plt.savefig(Path(path_cell_filtering_s, 'umap_doublet.png'), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_cell_filtering_s, 'umap_doublet.pdf'), bbox_inches="tight")

        with plt.rc_context():
            sc.pl.umap(
                adata_s,
                color=["log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
                wspace=0.3,
                ncols=2,
                show=False
            )
            plt.savefig(Path(path_cell_filtering_s, 'umap_total_counts_genes_mt.png'), bbox_inches="tight")    
            if args.pdf:
                plt.savefig(Path(path_cell_filtering_s, 'umap_total_counts_genes_mt.pdf'), bbox_inches="tight")    


    # add group column in adata.obs
    # if hasattr(samplesheet, 'group'):
    #     sample2group = dict(zip(samplesheet['sample'], samplesheet['group']))
    #     adata.obs['group'] = [sample2group.get(x, x) for x in adata.obs['sample']]

    # merge samples with the name in column 'merge' if exists
    if 'merge' in samplesheet.columns:
        if sum(samplesheet['merge'].notna())>0:
            ss1=samplesheet[samplesheet['merge'].notna()]
            sample2merge = dict(zip(ss1['sample'], ss1['merge']))
            adata.obs['sample'] = [sample2merge.get(x, x) for x in adata.obs['sample']]

    # save a filtered and normalized concated h5ad file
    adata.write_h5ad(Path(path_quant_qc, f'adata_filtered_normalized.h5ad'))

    # save analysis parameters into a json file
    with open(Path(path_quant_qc, 'parameters.json'), 'w') as file:
        params = {}
        params.update({"--h5ad": str(args.h5ad)})        
        params.update({"--samplesheet": str(args.samplesheet)})        
        params.update({"--min_genes": args.min_genes})        
        params.update({"--min_counts": args.min_counts})        
        if sum(args.max_genes) > 0: params.update({"--max_genes": args.max_genes})
        if sum(args.max_counts) > 0: params.update({"--max_counts": args.max_counts})
        if sum(args.min_gcounts) > 0: params.update({"--min_gcounts": args.min_gcounts})       
        params.update({"--min_cells": args.min_cells})        
        params.update({"--pct_mt": args.pct_mt})        
        if args.quantile_upper < 1: params.update({"--quantile_upper": args.quantile_upper})        
        if args.quantile_lower > 0: params.update({"--quantile_lower": args.quantile_lower})        
        params.update({"--iqr_coef": args.iqr_coef})        
        params.update({"--mt": args.mt})
        if args.find_doublets: 
            params.update({"--doublet_rate": args.doublet_rate})
        else:
            params.update({"--find_doublets": args.find_doublets})
        json.dump(params, file, indent=4)



if __name__ == "__main__":
    sys.exit(main())

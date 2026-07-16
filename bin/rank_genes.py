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
from scipy import sparse
import seaborn as sns
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
        help="If Specify a group name, compare with respect to this group.",
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
        help="Specify a column used to define cell-types for DEA between groups.",
    )
    parser.add_argument(
        "--celltypes",
        default=None,
        help="Specify a list cell-types for DEA between groups, e.g. 'celltype1,celltype2'.",
    )
    parser.add_argument(
        "--combine",
        help="Whether to combine all samples for marker gene identification.",
        action='store_true',
    )
    parser.add_argument(
        "--deseq2",
        help="Apply PyDESeq2 for pseudobulk differential expression analysis.",
        action='store_true',
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
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()

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
        if 'group' in adata.obs.columns:
            batch = 'group'
        elif 'plate' in adata.obs.columns:
            batch = 'plate'  
    else:
        batch = args.meta

    
    # differential expression analysis
    if args.deseq2: # Apply PyDESeq2 for pseudobulk differential expression analysis
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        if groups is None:
            groups = list(adata.obs[groupby].unique())
        groups.remove(args.reference)
        
        groups_to_test = [args.reference] + groups
        adata_sub = adata[adata.obs[groupby].isin(groups_to_test)].copy()
        adata_sub.obs[groupby] = pd.Categorical(adata_sub.obs[groupby], categories=groups_to_test)

        # Check whether each sample belongs to only one group
        sample_group_n = adata_sub.obs.groupby("sample")[groupby].nunique()
        if (sample_group_n > 1).any():
            bad_samples = sample_group_n[sample_group_n > 1].index.tolist()
            logger.error(f"There are samples belonging to multiple groups: {bad_samples}")
            sys.exit(2)

        if "counts" in adata_sub.layers:
            X = adata_sub.layers["counts"].copy()
        elif adata_sub.raw is not None:
            X = adata_sub.raw.X.copy()
        else:
            if 'log1p_total_counts' not in adata_sub.obs.columns:
                X = adata_sub.X.copy()
            else:
                logger.error(f"PyESeq2 requre raw counts for pseudobulk DEA, \
                             please provide raw counts in anndata.layers['counts'].")
                sys.exit(2)     
        if sparse.issparse(X):
            X = X.tocsr()
        else:
            X = sparse.csr_matrix(X) 

        # Normalise and log-transform for visualisation
        if 'log1p_total_counts' not in adata_sub.obs.columns:
            sc.pp.normalize_total(adata_sub, target_sum=1e4)
            sc.pp.log1p(adata_sub)                      

        # Create pseudobulk ID
        adata_sub.obs["pseudobulk_id"] = adata_sub.obs["sample"].astype(str)
        pseudobulk_id = pd.Categorical(adata_sub.obs["pseudobulk_id"])
        codes = pseudobulk_id.codes
        pb_names = pseudobulk_id.categories.astype(str)
        # Build pseudobulk aggregation matrix
        indicator = sparse.csr_matrix(
            (
                np.ones(adata_sub.n_obs),
                (codes, np.arange(adata_sub.n_obs))
            ),
            shape=(len(pb_names), adata_sub.n_obs)
        )
        # Aggregate cell counts into sample-level pseudobulk counts
        pb_counts = indicator @ X
        counts_df = pd.DataFrame(
            pb_counts.toarray(),
            index=pb_names,
            columns=adata_sub.var_names.astype(str)
        )
        counts_df = counts_df.round().astype(int)
        # create metadata
        metadata = (
            adata_sub.obs[["pseudobulk_id", "sample", groupby]]
            .drop_duplicates()
            .set_index("pseudobulk_id")
            .loc[counts_df.index]
        )
        metadata[groupby] = pd.Categorical(metadata[groupby], categories=groups_to_test)
        # filter genes
        keep_genes = counts_df.sum(axis=0) >= 10
        counts_df = counts_df.loc[:, keep_genes]

        # Run pyDESeq2
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata,
            design=f"~ {groupby}",
            refit_cooks=True,
            quiet=True,
        )
        dds.deseq2()

        # statistics and plots for each group vs reference
        for group in groups:
            stat_group_vs_ref = DeseqStats(
                dds,
                contrast=[groupby, group, args.reference],
                alpha=0.05,
            )
            stat_group_vs_ref.summary()

            res_group_vs_ref = stat_group_vs_ref.results_df.copy()
            res_group_vs_ref["gene"] = res_group_vs_ref.index
            # res_group_vs_ref["comparison"] = f"{group}_vs_{args.reference}"
            res_group_vs_ref = res_group_vs_ref.sort_values("padj")

            res_group_vs_ref.to_csv(Path(path_analysis, f"pydeseq2_{group}_vs_{args.reference}.csv"), index=False)

            # Ranking plot -----------------------------------------------------------
            # Remove genes with missing p-values
            res = res_group_vs_ref.dropna(subset=["pvalue", "log2FoldChange"]).copy()
            # Add -log10 p-value
            res["minus_log10_pvalue"] = -np.log10(res["pvalue"].clip(lower=1e-300))
            # Add significance category
            res["DE_status"] = "Not significant"
            res.loc[(res["padj"] < 0.05) & (res["log2FoldChange"] > 1), "DE_status"] = f"Up in {group}"
            res.loc[(res["padj"] < 0.05) & (res["log2FoldChange"] < -1), "DE_status"] = f"Up in {args.reference}"
            tops = res.sort_values("pvalue").head(args.n_genes).copy()
            tops = tops.sort_values("log2FoldChange", ascending=False)
            topgenes = tops["gene"].tolist()
            plt.figure(figsize=(9, 5))
            plt.bar(
                tops["gene"],
                tops["log2FoldChange"],
                color=["#4C72B0" if x > 0 else "#DD8452" for x in tops["log2FoldChange"]]
            )
            plt.axhline(0, color="black", linewidth=0.8)
            plt.xlabel("Gene")
            plt.ylabel("log2 fold change")
            plt.title(f"Top DE genes: {group} vs {args.reference}")
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()
            plt.savefig(Path(path_analysis, f"plot_genes_{group}_vs_{args.reference}.png"), bbox_inches="tight")
            if args.pdf:        
                plt.savefig(Path(path_analysis, f"plot_genes_{group}_vs_{args.reference}.pdf"), bbox_inches="tight")            
 
            # Dotplot for top DEGs -----------------------------------------------------
            groups_to_plot = [args.reference, group]
            plt.rcParams.update({
                "legend.fontsize": 'small',
            })
            adata_plot = adata_sub[adata_sub.obs[groupby].isin(groups_to_plot), topgenes].copy()
            adata_plot.obs[groupby] = pd.Categorical(adata_plot.obs[groupby], categories=groups_to_plot)
            sc.pl.dotplot(
                adata_plot,
                var_names=topgenes,
                groupby=groupby,
                standard_scale="var",
                dendrogram=False,
                swap_axes=True,
                show=False,
                figsize=(3, max(4, len(topgenes) * 0.3)), # Adjust height based on number of genes
            )
            plt.savefig(Path(path_analysis, f"dotplot_genes_{group}_vs_{args.reference}.png"), bbox_inches="tight")
            if args.pdf:
                plt.savefig(Path(path_analysis, f"dotplot_genes_{group}_vs_{args.reference}.pdf"), bbox_inches="tight")

            # Volcano plot --------------------------------------------------------------
            fig, ax = plt.subplots(figsize=(8, 8))
            sns.scatterplot(
                data=res,
                x="log2FoldChange",
                y="minus_log10_pvalue",
                hue="DE_status",
                s=40,
                edgecolor=None,
                ax=ax
            )
            ax.axvline(1, linestyle="--", color="grey", linewidth=0.8)
            ax.axvline(-1, linestyle="--", color="grey", linewidth=0.8)
            ax.axhline(-np.log10(0.05), linestyle="--", color="grey", linewidth=0.8)
            ax.set_xlabel("log2 fold change")
            ax.set_ylabel("-log10 p-value")
            ax.set_title(f"Volcano plot: {group} vs {args.reference}")
            ax.set_box_aspect(1) # Keep only the plotting panel square
            # Add some margin so gene labels are not clipped
            x_min = res["log2FoldChange"].min()
            x_max = res["log2FoldChange"].max()
            y_min = 0
            y_max = res["minus_log10_pvalue"].max()
            x_pad = (x_max - x_min) * 0.15
            y_pad = max(y_max * 0.15, 0.5)
            ax.set_xlim(x_min - x_pad, x_max + x_pad)
            ax.set_ylim(y_min, y_max + y_pad)
            # Label top 10 genes by p-value
            top_label = res.sort_values("pvalue").head(10)
            for _, row in top_label.iterrows():
                ax.text(
                    row["log2FoldChange"],
                    row["minus_log10_pvalue"],
                    row["gene"],
                    fontsize=8,
                    ha="left",
                    va="bottom",
                    clip_on=False
                )
            # Put legend outside the plot panel
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(
                handles=handles,
                labels=labels,
                title="DE status",
                loc="upper left",
                bbox_to_anchor=(1.03, 1),
                borderaxespad=0,
                frameon=False
            )
            # Reserve right-side space for legend
            fig.subplots_adjust(right=0.72)
            fig.savefig(Path(path_analysis, f"volcano_{group}_vs_{args.reference}.png"), bbox_inches="tight")
            if args.pdf:           
                fig.savefig(Path(path_analysis, f"volcano_{group}_vs_{args.reference}.pdf"), bbox_inches="tight")    

    elif groupby == 'group': # between conditions using Scanpy
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

    elif args.combine: #  one cluster vs rest for combined sample
        sc.tl.rank_genes_groups(
            adata, 
            groupby, 
            method=args.method, 
            groups=groups if groups else 'all', 
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
            plt.savefig(Path(path_analysis, f"plot_genes_{groupby}.png"), bbox_inches="tight")
            if args.pdf:
                plt.savefig(Path(path_analysis, f"plot_genes_{groupby}.pdf"), bbox_inches="tight")

        with plt.rc_context():
            sc.pl.rank_genes_groups_dotplot(
                adata, 
                n_genes=args.n_genes, 
                groups=groups,
            )
            plt.savefig(Path(path_analysis, f"dotplot_genes_{groupby}.png"), bbox_inches="tight")
            if args.pdf:
                plt.savefig(Path(path_analysis, f"dotplot_genes_{groupby}.pdf"), bbox_inches="tight")

        for gid in sorted(groups if groups else adata.obs[groupby].unique()):
            sc.get.rank_genes_groups_df(adata, group=gid).to_csv(
                Path(path_analysis, f'dea_{groupby}_{gid}_vs_{args.reference}.csv'), 
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
        if args.combine: params.update({"--combine": ''})
        if args.deseq2: params.update({"--deseq2": ''})
        json.dump(params, file, indent=4)



if __name__ == "__main__":
    sys.exit(main())

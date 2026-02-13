#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import argparse
import sys, json
from pathlib import Path
import util
import torch
import scvi


logger = util.get_named_logger('DEA_SCVI')



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
    parser.add_argument(
        "--covar_cat",
        type=util.stringlist,
        help="List of categorical covariates.",
        default=[],
    )
    parser.add_argument(
        "--covar_con",
        type=util.stringlist,
        help="List of continuous covariates.",
        default=[],
    )
    parser.add_argument(
        "--epochs",
        type=int,
        help="Set number of epochs for scvi training.",
        default=None,
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        help="Set batch size for scvi training.",
        default=None,
    )
    parser.add_argument(
        "--devices",
        type=int,
        help="Set number of cpus/gpus for scvi training.",
        default=None,
    )
    parser.add_argument(
        "--groupby",
        help="Group the cells by a category column.",
        default=None,
    )
    parser.add_argument(
        "--group1",
        type=util.stringlist,
        help="first groups of cells for DEA.",
        default=None,
    )
    parser.add_argument(
        "--group2",
        type=util.stringlist,
        help="second groups of cells for DEA.",
        default=None,
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
    # parser.add_argument(
    #     "--pairing",
    #     choices=['one-to-one', 'many-to-many'],
    #     help="how to pair the groups between group1 and group2",
    #     default=None,
    # )
    parser.add_argument(
        "--n_markers",
        type=int,
        help="Set number of top marker genes.",
        default=3,
    )
    parser.add_argument(
        "--deg_lfc",
        type=float,
        help="Set threshold of Log Fold Change for DEGs.",
        default=0,
    )
    parser.add_argument(
        "--deg_bayes",
        type=float,
        help="Set threshold of bayes factor for DEGs.",
        default=0,
    )
    parser.add_argument(
        "--deg_nzerosprop",
        type=float,
        help="Set threshold of proportion of non-zero expression cells in group1 for DEGs.",
        default=0,
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
    path_scvi_model = Path(path_analysis, "scvi_model")
    util.check_and_create_folder(path_scvi_model)

    adata = sc.read_h5ad(args.h5ad)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    elif adata.raw is not None:
        adata = adata.raw.to_adata()
    adata.layers["counts"] = adata.X.copy().tocsr()

    # Run Poisson-based HVG selection
    scvi.data.poisson_gene_selection(adata, layer='counts')
    adata = adata[:, adata.var["highly_variable"]].copy()

    batch_key = "plate" if "plate" in adata.obs.columns else "sample" # correct on plates for smart-seq data

    torch.set_float32_matmul_precision("high")

    covariate_kwargs = {}
    if args.covar_cat:
        covariate_kwargs["categorical_covariate_keys"] = args.covar_cat
    if args.covar_con:
        covariate_kwargs["continuous_covariate_keys"] = args.covar_con
    # adata.obs[batch_key] = adata.obs[batch_key].astype(str).astype('category')
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key=batch_key,
        **covariate_kwargs
    )
    model = scvi.model.SCVI(adata, gene_likelihood="nb")

    train_kwargs = {}
    if args.epochs is not None:
        train_kwargs["max_epochs"] = args.epochs
    if args.batch_size is not None:
        train_kwargs["batch_size"] = args.batch_size
    # if args.devices is not None:
    #     train_kwargs["devices"] = args.devices
    #     train_kwargs["strategy"] = "ddp_find_unused_parameters_true"
    #     ddp_active = True
    # else:
    #     ddp_active = False
    model.train(
        **train_kwargs,
        early_stopping=True,
        early_stopping_patience=20,
        early_stopping_monitor="elbo_validation"
    )

    model.save(path_scvi_model, overwrite=True)
    adata.obsm['X_scvi'] = model.get_latent_representation()
    adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size=1e4)
    # sc.pp.neighbors(adata, use_rep='X_scvi')
    # sc.tl.umap(adata)
 

    if args.meta == 'auto':
        # batch = 'group' if hasattr(adata.obs, 'group') else 'sample'
        batch = 'sample'
        if 'group' in adata.obs.columns:
            batch = 'group'
        elif 'plate' in adata.obs.columns:
            batch = 'plate' 
    else:
        batch = args.meta
    


    # with plt.rc_context():
    #     sc.pl.umap(
    #         adata,
    #         color=[args.groupby],
    #         # legend_loc="on data",
    #         show=False
    #     )
    #     plt.savefig(Path(path_analysis, f"umap_{args.groupby}.png"), bbox_inches="tight")
    #     if args.pdf:
    #         plt.savefig(Path(path_analysis, f"umap_{args.groupby}.pdf"), bbox_inches="tight")


    if args.groupby == 'group': # DEA between conditions

        if args.celltype_col: # DEA between conditions for each celltype
            celltypes = args.celltypes.split(',') if args.celltypes else sorted(adata.obs[args.celltype_col].unique()) 
            for celltype in celltypes:        
                adata_s = adata[adata.obs[args.celltype_col]==celltype].copy()   
                path_analysis_s = Path(path_analysis, f"celltype_{celltype}".replace(' ', '_').replace('/', '_'))
                util.check_and_create_folder(path_analysis_s)

                model_s = scvi.model.SCVI.load_query_data(adata_s, model)
                de_df_list = []
                for group1 in args.group1 if args.group1 else [None]:
                    for group2 in args.group2 if args.group2 else [None]:
                        de_df_c = model_s.differential_expression(
                            groupby=args.groupby, 
                            group1=group1, 
                            group2=group2,
                        )
                        de_df_list += [de_df_c]
                de_df = pd.concat(de_df_list, axis=0, ignore_index=False)
                eps = 1e-8
                de_df["lfc_log2"] = np.log2((de_df["raw_normalized_mean1"] + eps) / (de_df["raw_normalized_mean2"] + eps))

                de_df.to_csv(
                    Path(path_analysis_s, f'dea_csvi_{celltype}.csv'), 
                    index=True,
                )

                # marker genes doplot
                markers = {}
                for comp in de_df.comparison.unique():
                    comp_de_df = de_df.loc[de_df.comparison == comp]
                    comp_de_df = comp_de_df[comp_de_df["lfc_log2"] > args.deg_lfc]
                    comp_de_df = comp_de_df[comp_de_df["bayes_factor"] > args.deg_bayes]
                    comp_de_df = comp_de_df[comp_de_df["non_zeros_proportion1"] > args.deg_nzerosprop]
                    markers[comp] = comp_de_df.index.tolist()[:args.n_markers]

                sc.tl.dendrogram(adata_s, groupby=args.groupby, use_rep="X_scvi")
                markers = {k: v for k, v in markers.items() if len(v) > 0}
                with plt.rc_context():
                    sc.pl.dotplot(
                        adata_s,
                        markers,
                        groupby=args.groupby,
                        dendrogram=True,
                        color_map="Blues",
                        swap_axes=True,
                        use_raw=False,  
                        layer="counts",
                        standard_scale="var",
                    )
                    plt.savefig(Path(path_analysis_s, f"dotplot_{celltype}_{args.groupby}.png"), bbox_inches="tight")
                    if args.pdf:
                        plt.savefig(Path(path_analysis_s, f"dotplot_{celltype}_{args.groupby}.pdf"), bbox_inches="tight")

                # heatmap plot
                Ng = len(markers)
                with plt.rc_context():
                    sc.pl.heatmap(
                        adata_s,
                        markers,
                        groupby=args.groupby,
                        layer="scvi_normalized",
                        standard_scale="var",
                        dendrogram=True,
                        figsize=((Ng+3)*0.8, adata_s.n_obs*0.004),
                    )  
                    plt.savefig(Path(path_analysis_s, f"heatmap_{celltype}_{args.groupby}.png"), bbox_inches="tight")
                    if args.pdf:
                        plt.savefig(Path(path_analysis_s, f"heatmap_{celltype}_{args.groupby}.pdf"), bbox_inches="tight")     

        else: # DEA between conditions for all cells
            de_df_list = []
            for group1 in args.group1 if args.group1 else [None]:
                for group2 in args.group2 if args.group2 else [None]:
                    de_df_c = model.differential_expression(
                        groupby=args.groupby, 
                        group1=group1, 
                        group2=group2,
                    )
                    de_df_list += [de_df_c]
            de_df = pd.concat(de_df_list, axis=0, ignore_index=False)
            eps = 1e-8
            de_df["lfc_log2"] = np.log2((de_df["raw_normalized_mean1"] + eps) / (de_df["raw_normalized_mean2"] + eps))        

            # for comp in de_df['comparison'].unique():
            #     de_df[de_df['comparison']==comp].to_csv(
            #         Path(path_analysis, f'dea_{comp}.csv'), 
            #         index=True,
            #     )

            # write all dea comparisons into one csv file
            de_df.to_csv(
                Path(path_analysis, f'dea_csvi.csv'), 
                index=True,
            )

            # marker genes doplot
            markers = {}
            for comp in de_df.comparison.unique():
                comp_de_df = de_df.loc[de_df.comparison == comp]
                comp_de_df = comp_de_df[comp_de_df["lfc_log2"] > args.deg_lfc]
                comp_de_df = comp_de_df[comp_de_df["bayes_factor"] > args.deg_bayes]
                comp_de_df = comp_de_df[comp_de_df["non_zeros_proportion1"] > args.deg_nzerosprop]
                markers[comp] = comp_de_df.index.tolist()[:args.n_markers]

            sc.tl.dendrogram(adata, groupby=args.groupby, use_rep="X_scvi")
            markers = {k: v for k, v in markers.items() if len(v) > 0}
            with plt.rc_context():
                sc.pl.dotplot(
                    adata,
                    markers,
                    groupby=args.groupby,
                    dendrogram=True,
                    color_map="Blues",
                    swap_axes=True,
                    use_raw=False,  
                    layer="counts",
                    standard_scale="var",
                )
                plt.savefig(Path(path_analysis, f"dotplot_{args.groupby}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_analysis, f"dotplot_{args.groupby}.pdf"), bbox_inches="tight")

            # heatmap plot
            Ng = len(markers)
            with plt.rc_context():
                sc.pl.heatmap(
                    adata,
                    markers,
                    groupby=args.groupby,
                    layer="scvi_normalized",
                    standard_scale="var",
                    dendrogram=True,
                    figsize=((Ng+3)*0.8, adata.n_obs*0.004),
                )  
                plt.savefig(Path(path_analysis, f"heatmap_{args.groupby}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_analysis, f"heatmap_{args.groupby}.pdf"), bbox_inches="tight")  

    else: # grp1 vs grp2 for each sample/group
        for sid in sorted(adata.obs[batch].unique()):
            adata_s = adata[adata.obs[batch]==sid].copy()   
            path_analysis_s = Path(path_analysis, f"{batch}_{sid}")
            util.check_and_create_folder(path_analysis_s)

            model_s = scvi.model.SCVI.load_query_data(adata_s, model)
            de_df_list = []
            for group1 in args.group1 if args.group1 else [None]:
                for group2 in args.group2 if args.group2 else [None]:
                    de_df_c = model_s.differential_expression(
                        groupby=args.groupby, 
                        group1=group1, 
                        group2=group2,
                    )
                    de_df_list += [de_df_c]
            de_df = pd.concat(de_df_list, axis=0, ignore_index=False)
            eps = 1e-8
            de_df["lfc_log2"] = np.log2((de_df["raw_normalized_mean1"] + eps) / (de_df["raw_normalized_mean2"] + eps))

            de_df.to_csv(
                Path(path_analysis_s, f'dea_csvi_{sid}.csv'), 
                index=True,
            )

            # marker genes doplot
            markers = {}
            for comp in de_df.comparison.unique():
                comp_de_df = de_df.loc[de_df.comparison == comp]
                comp_de_df = comp_de_df[comp_de_df["lfc_log2"] > args.deg_lfc]
                comp_de_df = comp_de_df[comp_de_df["bayes_factor"] > args.deg_bayes]
                comp_de_df = comp_de_df[comp_de_df["non_zeros_proportion1"] > args.deg_nzerosprop]
                markers[comp] = comp_de_df.index.tolist()[:args.n_markers]

            sc.tl.dendrogram(adata_s, groupby=args.groupby, use_rep="X_scvi")
            markers = {k: v for k, v in markers.items() if len(v) > 0}
            with plt.rc_context():
                sc.pl.dotplot(
                    adata_s,
                    markers,
                    groupby=args.groupby,
                    dendrogram=True,
                    color_map="Blues",
                    swap_axes=True,
                    use_raw=False,  
                    layer="counts",
                    standard_scale="var",
                )
                plt.savefig(Path(path_analysis_s, f"dotplot_{sid}_{args.groupby}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_analysis_s, f"dotplot_{sid}_{args.groupby}.pdf"), bbox_inches="tight")

            # heatmap plot
            Ng = len(markers)
            with plt.rc_context():
                sc.pl.heatmap(
                    adata_s,
                    markers,
                    groupby=args.groupby,
                    layer="scvi_normalized",
                    standard_scale="var",
                    dendrogram=True,
                    figsize=((Ng+3)*0.8, adata_s.n_obs*0.004),
                )  
                plt.savefig(Path(path_analysis_s, f"heatmap_{sid}_{args.groupby}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_analysis_s, f"heatmap_{sid}_{args.groupby}.pdf"), bbox_inches="tight")            


    # save analysis parameters into a json file
    with open(Path(path_analysis, 'parameters.json'), 'w') as file:
        params = {}
        params.update({"--h5ad": str(args.h5ad)})      
        params.update({"--groupby": args.groupby})
        params.update({"--group1": args.group1})
        params.update({"--group2": args.group2})
        params.update({"--n_markers": args.n_markers})
        params.update({"--meta": args.meta})        
        if args.covar_cat: params.update({"--covar_cat": args.covar_cat})
        if args.covar_con: params.update({"--covar_con": args.covar_con})
        if args.epochs: params.update({"--epochs": args.epochs})
        if args.batch_size: params.update({"--batch_size": args.batch_size})
        if args.devices: params.update({"--devices": args.devices})
        json.dump(params, file, indent=4)



if __name__ == "__main__":
    sys.exit(main())

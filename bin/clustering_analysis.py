#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
# import torch.distributed as dist
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'
# IS_MAIN_PROCESS = os.environ.get("GLOBAL_RANK", "0") == "0"

import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt
import argparse
import sys, json
from pathlib import Path
import numpy as np
import util
# import torch
# import scvi


logger = util.get_named_logger('CLUSTERING')



# import torch.distributed as dist

# def is_main_process():
#     # Torch DDP sets this automatically
#     if dist.is_available() and dist.is_initialized():
#         return dist.get_rank() == 0

#     # Fall back to environment vars if dist not yet initialized
#     for key in ["RANK", "SLURM_PROCID", "LOCAL_RANK", "GLOBAL_RANK"]:
#         if key in os.environ:
#             return os.environ[key] in ("0", 0)
#     return True

# IS_MAIN_PROCESS = is_main_process()


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
        "--normalize",
        help="Whether to apply normalization to the data.",
        action='store_true',
        # action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--keep_doublets",
        help="Whether to filter out the cells called as doublets.",
        action='store_true',
    )       
    parser.add_argument(
        "--regress",
        help="Whether to regress out the variations from the total counts and the percentage of mitochondrial genes expressed.",
        action='store_true',
    )
    parser.add_argument(
        "--scale",
        help="Whether to scale the expression to have zero mean and unit variance.",
        action='store_true',
    )        
    parser.add_argument(
        "--n_comps",
        type=int,
        help="Set the number of components for PCA dimensionality reduction.",
        default=50,
    )
    parser.add_argument(
        "--n_neighbors",
        type=int,
        help="Set the number of neighbors for nearest neighbor graph constuction.",
        default=15,
    )     
    parser.add_argument(
        "--resolutions",
        type=util.floatlist,
        help="Resolution is used to control number of clusters.",
        default=[0.02, 0.05, 0.1, 0.5],
    )
    parser.add_argument(
        "--integrate",
        default=None,
        choices=['bbknn', 'harmony', 'scanorama', 'scvi'],
        help="Choose a method for data integration",
    )
    parser.add_argument(
        "--min_cluster_size",
        type=int,
        help="Set minimal cluster size for filtering out small clusters.",
        default=0,
    )
    parser.add_argument(
        "--min_cluster_pct",
        type=float,
        help="Set minimal cluster percentage number for filtering out small clusters.",
        default=0,
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
    path_clustering = Path(args.outdir)
    util.check_and_create_folder(path_clustering)

    adata = sc.read_h5ad(args.h5ad)

    # save the original counts
    if "counts" not in adata.layers: 
        adata.layers["counts"] = adata.X.copy()
        
    # Normalization
    if args.normalize:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    if not adata.uns.get('log1p'): # to fix issue in scanpy function
        adata.uns['log1p'] = {'base': None}

    # remove doublets before clustering
    if not args.keep_doublets:
        if hasattr(adata.obs, 'predicted_doublet'):
            adata = adata[~adata.obs['predicted_doublet']]

    batch_key = 'plate' if hasattr(adata.obs, 'plate') else 'sample' # correct on plates for smart-seq data

    # Feature selection and dimensionality reduction
    sc.pp.highly_variable_genes(
        adata,
        flavor='seurat', 
        min_mean=0.0125, 
        max_mean=1000000000, 
        min_disp=0.5, 
        max_disp=50, 
        n_bins=20, 
        batch_key=batch_key
    )
    #sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")

    # save the raw counts
    #s adata.raw = adata
    #s adata = adata[:, adata.var.highly_variable]


    if args.integrate == 'scvi':
        import torch # only if cpu support AVX instructions
        import scvi
        path_scvi_model = Path(path_clustering, "scvi_model")
        util.check_and_create_folder(path_scvi_model)
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
        model = scvi.model.SCVI(adata)

        train_kwargs = {}
        if args.epochs is not None:
            train_kwargs["max_epochs"] = args.epochs
        if args.batch_size is not None:
            train_kwargs["batch_size"] = args.batch_size
        # if args.devices is not None:
        #     train_kwargs["devices"] = args.devices
        #     train_kwargs["strategy"] = "ddp_find_unused_parameters_true"
        model.train(**train_kwargs)
        
        # if args.devices is not None:
        # if not IS_MAIN_PROCESS: sys.exit(0)
        # if IS_MAIN_PROCESS:
        model.save(path_scvi_model, overwrite=True)
        adata.obsm['X_scvi'] = model.get_latent_representation()
        adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size=1e4)
        sc.pp.neighbors(adata, use_rep='X_scvi')
        sc.tl.umap(adata)

    else:
        # Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed
        if args.regress:
            if not(hasattr(adata.obs, 'total_counts') and hasattr(adata.obs, 'pct_counts_mt')):
                adata.var['mt'] = adata.var_names.str.startswith('MT-')
                sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)
                # sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)    
            sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

        # scale the expression to have zero mean and unit variance
        if args.scale:
            sc.pp.scale(adata, max_value=10)

        # Dimensionality reduction
        n_comps = min((min(adata.X[:, adata.var["highly_variable"].values].shape)-1), 50)
        sc.tl.pca(
            adata, 
            n_comps=n_comps
            # chunked=False,
            # zero_center=False, 
            # svd_solver='arpack'
        )

        # perform data integration
        n_pcs = adata.obsm['X_pca'].shape[1]
        if args.integrate == 'bbknn':
            sc.external.pp.bbknn(adata, batch_key=batch_key, n_pcs=n_pcs)
        elif args.integrate == 'harmony':
            mask = ~np.isnan(adata.obsm["X_pca"]).any(axis=1)
            if not np.all(mask):
                logger.warning(f"Removing {np.sum(~mask)} cells with NaN PCA embeddings before Harmony.")
                adata = adata[mask].copy()
            sc.external.pp.harmony_integrate(adata, batch_key)
        elif args.integrate == 'scanorama':
            sc.external.pp.scanorama_integrate(adata, batch_key)
            
        # find nearest neighbor graph constuction
        pca_rep = 'X_pca'
        if args.integrate == 'harmony':
            pca_rep = 'X_pca_harmony'
        elif args.integrate == 'scanorama':
            pca_rep = 'X_scanorama'
        if args.integrate != 'bbknn':
            sc.pp.neighbors(
                adata, 
            n_neighbors=args.n_neighbors, 
                n_pcs=n_pcs,
                knn=True, 
                method='umap', 
                metric='euclidean',
                use_rep=pca_rep
            )
        sc.tl.umap(adata)


    # perform clustering using Leiden graph-clustering method
    for res in args.resolutions:
        sc.tl.leiden(
            adata, n_iterations=2, 
            key_added=f"leiden_res_{res:4.2f}", resolution=res
        )

    # Filter Out Small Clusters
    if args.min_cluster_size > 0:
        res_columns = [f"leiden_res_{res:4.2f}" for res in args.resolutions]
        keep_mask = np.array([True] * adata.n_obs)
        for col in res_columns:
            if col in adata.obs.columns:
                counts = adata.obs[col].value_counts()
                small_clusters = counts[counts < args.min_cluster_size].index.tolist()
                is_in_small_cluster = adata.obs[col].isin(small_clusters)
                keep_mask = keep_mask & (~is_in_small_cluster)
        adata = adata[keep_mask].copy()

    if args.min_cluster_pct > 0:
        res_columns = [f"leiden_res_{res:4.2f}" for res in args.resolutions]
        keep_mask = np.array([True] * adata.n_obs)
        for col in res_columns:
            if col in adata.obs.columns:
                pcts = adata.obs[col].value_counts(normalize=True) * 100
                small_clusters = pcts[pcts < args.min_cluster_pct].index.tolist()
                is_in_small_cluster = adata.obs[col].isin(small_clusters)
                keep_mask = keep_mask & (~is_in_small_cluster)
        adata = adata[keep_mask].copy()

    # save the AnnData into a h5ad file 
    adata.write_h5ad(Path(path_clustering, 'adata_clustering.h5ad'))

    if args.meta == 'auto':
        # batch = 'group' if hasattr(adata.obs, 'group') else 'sample'
        batch = 'sample'
        if hasattr(adata.obs, 'group'):
            batch = 'group'
        elif hasattr(adata.obs, 'plate'):
            batch = 'plate' 
    else:
        batch = args.meta
    
    # sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=args.resolution)
    for sid in sorted(adata.obs[batch].unique()):
        adata_s = adata[adata.obs[batch]==sid]   
        path_clustering_s = Path(path_clustering, f"{batch}_{sid}")
        util.check_and_create_folder(path_clustering_s)
        for res in args.resolutions:
            with plt.rc_context():
                sc.pl.umap(
                    adata_s,
                    color=[f"leiden_res_{res:4.2f}"],
                    # legend_loc="on data",
                    show=False
                )
                plt.savefig(Path(path_clustering_s, f"umap_leiden_res_{res:4.2f}.png"), bbox_inches="tight")
                if args.pdf:
                    plt.savefig(Path(path_clustering_s, f"umap_leiden_res_{res:4.2f}.pdf"), bbox_inches="tight")

    
    # stacked proportion bar plot showing all samples for each resolution      
    for res in args.resolutions:
        sc.pl.umap(adata, color=f'leiden_res_{res:4.2f}', show=False)
        path_res = Path(path_clustering, f"resolution_{res:4.2f}")
        util.check_and_create_folder(path_res)
        n_cluster = len(adata.obs[f'leiden_res_{res:4.2f}'].unique())+1
        ncol = min((n_cluster//20 + min(n_cluster%20, 1)), 3)
        colors = list(adata.uns[f'leiden_res_{res:4.2f}_colors'])
        with plt.rc_context():
            prop = pd.crosstab(adata.obs[f'leiden_res_{res:4.2f}'],adata.obs[batch], normalize='columns').T.plot(kind='bar', stacked=True, color=colors)
            prop.legend(bbox_to_anchor=(1.05+(args.fontsize-10)/50+ncol*0.17, 1.02), loc='upper right', ncol=ncol)
            plt.savefig(Path(path_res, f"prop_leiden_res_{res:4.2f}.png"), bbox_inches="tight")
            if args.pdf:
                plt.savefig(Path(path_res, f"prop_leiden_res_{res:4.2f}.pdf"), bbox_inches="tight")


    # save analysis parameters into a json file
    with open(Path(path_clustering, 'parameters.json'), 'w') as file:
        params = {}
        params.update({"--h5ad": str(args.h5ad)})        
        params.update({"--n_comps": args.n_comps})        
        params.update({"--n_neighbors": args.n_neighbors})        
        params.update({"--resolutions": args.resolutions})        
        if args.integrate: params.update({"--integrate": args.integrate})        
        if args.min_cluster_size > 0: params.update({"--min_cluster_size": args.min_cluster_size})    
        if args.min_cluster_pct > 0: params.update({"--min_cluster_pct": args.min_cluster_pct})    
        params.update({"--meta": args.meta})        
        if args.normalize: params.update({"--normalize": args.normalize})
        if args.regress: params.update({"--regress": args.regress})
        if args.scale: params.update({"--scale": args.scale})
        if args.keep_doublets: params.update({"--keep_doublets": args.keep_doublets})
        if args.covar_cat: params.update({"--covar_cat": args.covar_cat})
        if args.covar_con: params.update({"--covar_con": args.covar_con})
        if args.epochs: params.update({"--epochs": args.epochs})
        if args.batch_size: params.update({"--batch_size": args.batch_size})
        # if args.devices: params.update({"--devices": args.devices})
        json.dump(params, file, indent=4)



if __name__ == "__main__":
    sys.exit(main())

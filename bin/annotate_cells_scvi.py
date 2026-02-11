#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from matplotlib import pyplot as plt
import argparse
import sys, json
from pathlib import Path
import util
import scvi
import torch

logger = util.get_named_logger('ANNOTATE_CELLS')

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
        help="Input query anndata data file.",
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
        "--h5ad_ref",
        metavar="REF_H5AD",
        type=Path,
        help="Input reference anndata data file.",
        # required=True,
    )
    parser.add_argument(
        "--model_path",
        metavar="MODELFILE",
        type=Path,
        help="Path to scANVI model files.",
    )
    # parser.add_argument(
    #     "--model_cls",
    #     default='scvi',
    #     choices=['scvi', 'scanvi'],
    #     help="The model class to be trained.",
    # )
    parser.add_argument(
        "--batch_key",
        default='sample',
        help="Specify a batch key for reference data.",
    )
    parser.add_argument(
        "--early_stop",
        help="Whether to stop training based on validation metrics.",
        action='store_true',
    )
    parser.add_argument(
        "--label_key",
        default=None,
        help="Specify a label key of cell types for reference data.",
    )
    parser.add_argument(
        "--n_top_genes",
        type=int,
        help="Set number of highly-variable genes to keep.",
        default=2000,
    )
    parser.add_argument(
        "--min_label_pct",
        type=float,
        help="Set minimal label percentage number for filtering out small labels.",
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
        "--scvi_epochs",
        type=int,
        help="Set number of epochs for scvi training.",
        default=None,
    )
    parser.add_argument(
        "--scanvi_epochs",
        type=int,
        help="Set number of epochs for scanvi training.",
        default=None,
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        help="Set batch size for scvi training.",
        default=None,
    )
    parser.add_argument(
        "--n_samples_pl",
        type=int,
        help="Set number of samples per label for scanvi training.",
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

    if args.model_path and not args.model_path.is_dir():
        logger.error(f"The given input file {args.model_path} was not found!")
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
    path_annotation = Path(args.outdir)
    util.check_and_create_folder(path_annotation)

    torch.set_float32_matmul_precision("high")

    adata = sc.read_h5ad(args.h5ad)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    elif adata.raw is not None:
        adata = adata.raw.to_adata()

    scanvi_ref = None
    if args.h5ad_ref:
        adata_ref = sc.read_h5ad(args.h5ad_ref)
        # Clean up Inf and NaN values in X
        if sp.issparse(adata_ref.X):
            adata_ref.X = adata_ref.X.toarray()
        adata_ref.X = np.nan_to_num(adata_ref.X, nan=0.0, posinf=0.0, neginf=0.0)
        # Normalize/log to avoid huge magnitudes producing inf later
        sc.pp.normalize_total(adata_ref, target_sum=1e4, inplace=True)
        sc.pp.log1p(adata_ref)

        sc.pp.highly_variable_genes(
            adata_ref, 
            n_top_genes=args.n_top_genes, 
            batch_key=args.batch_key,
            subset=True
        )
        # sc.pp.highly_variable_genes(
        #     adata_ref,
        #     flavor='seurat', 
        #     min_mean=0.0125, 
        #     max_mean=1000000000, 
        #     min_disp=0.5, 
        #     max_disp=50, 
        #     n_bins=20, 
        #     batch_key=args.batch_key,
        #     subset=True
        # )

        common_genes = adata_ref.var_names.intersection(adata.var_names)
        adata_ref = adata_ref[:, common_genes].copy()
        adata = adata[:, common_genes].copy()
        # adata = adata[:, adata_ref.var_names].copy()

        if "counts" not in adata_ref.layers: 
            adata_ref.layers["counts"] = adata_ref.X.copy()
        scvi.model.SCVI.setup_anndata(adata_ref, batch_key=args.batch_key, layer="counts")
        scvi_ref = scvi.model.SCVI(
            adata_ref,
            use_layer_norm="both",
            use_batch_norm="none",
            encode_covariates=True,
            dropout_rate=0.2,
            n_layers=2,
        )

        train_kwargs = {}
        if args.scvi_epochs is not None:
            train_kwargs["max_epochs"] = args.scvi_epochs
        if args.batch_size is not None:
            train_kwargs["batch_size"] = args.batch_size
        # if args.devices is not None:
        #     train_kwargs["devices"] = args.devices
        #     train_kwargs["strategy"] = "ddp_find_unused_parameters_true"
        if args.early_stop:
            train_kwargs["early_stopping"] = True
            train_kwargs["early_stopping_patience"] = 20
            train_kwargs["early_stopping_min_delta"] = 1e-3        
        scvi_ref.train(**train_kwargs)

        # if not IS_MAIN_PROCESS: sys.exit(0)

        # adata_ref.obsm['X_scvi'] = scvi_ref.get_latent_representation()
        scanvi_ref = scvi.model.SCANVI.from_scvi_model(
            scvi_ref,
            unlabeled_category='Unknown',
            labels_key=args.label_key,
        )

        train_kwargs = {}
        if args.scanvi_epochs is not None:
            train_kwargs["max_epochs"] = args.scanvi_epochs
        if args.batch_size is not None:
            train_kwargs["batch_size"] = args.batch_size
        if args.n_samples_pl is not None:
            train_kwargs["n_samples_per_label"] = args.n_samples_pl
        # if args.devices is not None:
        #     train_kwargs["devices"] = args.devices
        #     train_kwargs["strategy"] = "ddp_find_unused_parameters_true"
        if args.early_stop:
            train_kwargs["early_stopping"] = True
            train_kwargs["early_stopping_patience"] = 20
            train_kwargs["early_stopping_min_delta"] = 1e-3 
        scanvi_ref.train(**train_kwargs)  
        # scanvi_ref.train(max_epochs=20, n_samples_per_label=100, devices=args.devices)

        # if not IS_MAIN_PROCESS: sys.exit(0)

        # scanvi_ref.save(Path(path_annotation, 'scanvi_model'), overwrite=True)
    else:
        scanvi_ref = scvi.model.SCANVI.load(args.model_path)
        batch_key = scanvi_ref.registry.get('setup_args').get('batch_key')
        if batch_key:
            if batch_key not in adata.obs.columns:
                adata.obs[batch_key] = adata.obs[args.batch_key].copy()
        else:
            logger.error(f"Can't find 'batch_key' in the model!")
            sys.exit(2)           

    # transfer learning on query
    scvi.model.SCANVI.prepare_query_anndata(adata, scanvi_ref)
    scanvi_query = scvi.model.SCANVI.load_query_data(adata, scanvi_ref)
    scanvi_query.train(
        max_epochs=100,
        plan_kwargs={"weight_decay": 0.0},
        check_val_every_n_epoch=10,
    )
    adata.obsm['X_scanvi'] = scanvi_query.get_latent_representation()
    adata.obs['scanvi_label'] = scanvi_query.predict()
    # adata.obs['scanvi_prob'] = scanvi_query.predict_proba()
    probs = scanvi_query.predict(adata=adata, soft=True)
    adata.obs['scanvi_prob'] = np.asarray(probs).max(axis=1)
    sc.pp.neighbors(adata, use_rep='X_scanvi')
    sc.tl.umap(adata)


    # Filter Out Small labels
    if args.min_label_pct > 0:
        keep_mask = np.array([True] * adata.n_obs)
        pcts = adata.obs['scanvi_label'].value_counts(normalize=True) * 100
        small_clusters = pcts[pcts < args.min_label_pct].index.tolist()
        is_in_small_cluster = adata.obs['scanvi_label'].isin(small_clusters)
        keep_mask = keep_mask & (~is_in_small_cluster)
        adata = adata[keep_mask].copy()

    # save the AnnData into a h5ad file
    adata.write_h5ad(Path(path_annotation, 'adata_annotation.h5ad'))

    # plot UMAPs to visualise the predicted cell-types
    if args.meta == 'auto':
        # batch = 'group' if hasattr(adata.obs, 'group') else 'sample'
        batch = 'sample'
        if hasattr(adata.obs, 'group'):
            batch = 'group'
        elif hasattr(adata.obs, 'plate'):
            batch = 'plate'     
    else:
        batch = args.meta

    label_type = 'scanvi_label'
    for sid in sorted(adata.obs[batch].unique()):
        adata_s = adata[adata.obs[batch]==sid]   
        path_annotation_s = Path(path_annotation, f"{batch}_{sid}")
        util.check_and_create_folder(path_annotation_s)
        n_cluster = len(adata_s.obs[label_type].unique())+1
        ncol = min((n_cluster//20 + min(n_cluster%20, 1)), 3)
        with plt.rc_context({"figure.dpi": 100}):
            # sc.pl.umap(
            #     adata_s,
            #     color=label_type,
            #     # legend_loc="on data",
            #     show=False
            # )

            ax = sc.pl.umap(adata_s, color=label_type, show=False)
            ax.set_aspect("equal", adjustable="box")
            fig = ax.figure
            leg = ax.legend_
            if leg is not None:
                bbox = leg.get_window_extent()
                width_px = bbox.width
                if width_px > 500:
                    leg.remove()
                    ax.legend(bbox_to_anchor=(0.5, -0.15), loc="upper center", ncol=min(ncol+1, 6))

            plt.savefig(Path(path_annotation_s, f"umap_cell_type.png"), bbox_inches="tight", dpi=150)    
            if args.pdf:
                plt.savefig(Path(path_annotation_s, f"umap_cell_type.pdf"), bbox_inches="tight")    

        with plt.rc_context():
            sc.pl.umap(
                adata_s,
                color='scanvi_prob',
                # legend_loc="on data",
                show=False
            )
            plt.savefig(Path(path_annotation_s, f"umap_conf_score.png"), bbox_inches="tight")          
            if args.pdf:
                plt.savefig(Path(path_annotation_s, f"umap_conf_score.pdf"), bbox_inches="tight")          


    # stacked proportion bar plot to compare between batches
    sc.pl.umap(adata, color=label_type, show=False) # get adata.uns['_colors']     
    n_cluster = len(adata.obs[label_type].unique())+1
    ncol = min((n_cluster//20 + min(n_cluster%20, 1)), 3)
    with plt.rc_context():
        ax = pd.crosstab(
            adata.obs[label_type], 
            adata.obs[batch], 
            normalize='columns').T.plot(kind='bar', stacked=True, color=adata.uns[f"{label_type}_colors"])
        # prop.legend(bbox_to_anchor=(1.4+(args.fontsize-10)/50+ncol*0.17, 1.02), loc='upper right', ncol=ncol)

        # leg = plt.legend(
        #     bbox_to_anchor=(1.4 + (args.fontsize - 10) / 50 + ncol * 0.17, 1.02),
        #     loc="upper right",
        #     ncol=ncol,
        # )
        leg = ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.05, 1.0),   # just outside the axes
            borderaxespad=0,
            ncol=ncol,
        )
        bbox = leg.get_window_extent()
        width = bbox.width
        if width > 500:
            leg.remove()
            plt.legend(bbox_to_anchor=(0.5, -0.4), loc="upper center", ncol=min(ncol+1, 6))

        plt.savefig(Path(path_annotation, f"prop_{label_type}.png"), bbox_inches="tight")    
        if args.pdf:
            plt.savefig(Path(path_annotation, f"prop_{label_type}.pdf"), bbox_inches="tight")    


    # save analysis parameters into a json file
    with open(Path(path_annotation, 'parameters.json'), 'w') as file:
        params = {}
        params.update({"--h5ad": str(args.h5ad)})        
        if args.h5ad_ref: params.update({"--h5ad_ref": str(args.h5ad_ref)})        
        if args.model_path: 
            params.update({"--model_path": str(args.model_path)})
        params.update({"--batch_key": str(args.batch_key)}) 
        params.update({"--label_key": str(args.label_key)}) 
        params.update({"--n_top_genes": str(args.n_top_genes)})
        if args.early_stop: params.update({"--early_stop": args.early_stop})  
        if args.min_label_pct > 0: params.update({"--min_label_pct": args.min_label_pct})  
        params.update({"--meta": args.meta})        
        if args.scvi_epochs: params.update({"--scvi_epochs": args.scvi_epochs})
        if args.scanvi_epochs: params.update({"--scanvi_epochs": args.scanvi_epochs})
        if args.batch_size: params.update({"--batch_size": args.batch_size})
        # if args.devices: params.update({"--devices": args.devices})
        if args.n_samples_pl: params.update({"--n_samples_pl": args.n_samples_pl})
        json.dump(params, file, indent=4)



if __name__ == "__main__":
    sys.exit(main())

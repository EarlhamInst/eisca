#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'
# os.mkdir('/.celltypist')
# os.chmod('/.celltypist', 0o755)
# os.environ['CELLTYPIST_PATH'] = '/.celltypist'
# os.environ['CELLTYPIST_PATH'] = os.path.expanduser('~/.celltypist')

import scanpy as sc
import celltypist as ct
import pandas as pd
from matplotlib import pyplot as plt
import argparse
import sys, json
from pathlib import Path
import util

logger = util.get_named_logger('ANNOTATE_CELLS')


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
        "--model_file",
        metavar="MODELFILE",
        type=Path,
        help="Path to CellTypist model file.",
    )
    parser.add_argument(
        "--model",
        default='Immune_All_Low.pkl',
        help="Specify a CellTypist model name, igored if a model file specified.",
    )    
    parser.add_argument(
        "--mode",
        default='best match',
        choices=['best match', 'prob match'],
        help="The cell-type prediction mode.",
    )
    parser.add_argument(
        "--p_thres",
        type=float,
        default=0.5,
        help="Probability threshold for the multi-label classification. Ignored if mode is 'best match'.",
    )              
    parser.add_argument(
        "--no_majority_voting",
        help="Whether to disable the majority voting classifier.",
        action='store_true',
    )
    parser.add_argument(
        "--update_models",
        help="Whether to update CellTypist models.",
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

    if args.model_file and not args.model_file.is_file():
        logger.error(f"The given input file {args.model_file} was not found!")
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

    adata = sc.read_h5ad(args.h5ad)

    # update CellTypist cell-type models
    if args.update_models:
        ct.models.models_description(on_the_fly = True)
    

    # perform cell-type annotation
    # model = ct.models.Model.load(model = 'Immune_All_Low.pkl')
    model = str(args.model_file) if args.model_file else args.model
    majority_voting = False if args.no_majority_voting else True
    predictions = ct.annotate(
        adata, 
        model=model, 
        majority_voting=majority_voting, 
        mode=args.mode,
        p_thres=args.p_thres,
    )
    adata = predictions.to_adata()


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

    label_type = 'majority_voting' if hasattr(adata.obs, 'majority_voting') else 'predicted_labels'
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
                    ax.legend(bbox_to_anchor=(0.5, -0.2), loc="upper center", ncol=min(ncol+1, 6))

            plt.savefig(Path(path_annotation_s, f"umap_cell_type.png"), bbox_inches="tight")    
            if args.pdf:
                plt.savefig(Path(path_annotation_s, f"umap_cell_type.pdf"), bbox_inches="tight")    

        with plt.rc_context():
            sc.pl.umap(
                adata_s,
                color='conf_score',
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
        # prop = pd.crosstab(adata.obs[label_type], adata.obs[batch], normalize='columns').T.plot(kind='bar', stacked=True, color=adata.uns[f"{label_type}_colors"])
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
            leg.remove()  # remove old one
            plt.legend(bbox_to_anchor=(0.5, -0.4), loc="upper center", ncol=min(ncol+1, 6))

        plt.savefig(Path(path_annotation, f"prop_{label_type}.png"), bbox_inches="tight")    
        if args.pdf:
            plt.savefig(Path(path_annotation, f"prop_{label_type}.pdf"), bbox_inches="tight")  


    # save analysis parameters into a json file
    with open(Path(path_annotation, 'parameters.json'), 'w') as file:
        params = {}
        params.update({"--h5ad": str(args.h5ad)})        
        if args.model_file: 
            params.update({"--model_file": str(args.model_file)})
        else:
            params.update({"--model": args.model})   
        params.update({"--mode": args.mode})        
        params.update({"--p_thres": args.p_thres})        
        params.update({"--meta": args.meta})        
        if args.no_majority_voting: params.update({"--no_majority_voting": args.no_majority_voting})
        if args.update_models: params.update({"--update_models": args.update_models})
        json.dump(params, file, indent=4)



if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import scanpy as sc
import celltypist as ct
import pandas as pd
from matplotlib import pyplot as plt
import argparse
import sys
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
        "--model_filename",
        default='celltypist_model.pkl',
        help="Specify a CellTypist model name.",
    )
    parser.add_argument(
        "--labels",
        default='cell_type',
        help="Specify a column of cell-type from cell metadata of Anndata.",
    )         
    parser.add_argument(
        "--l2c",
        type=float,
        default=1.0,
        help="Inverse of L2 regularization strength for traditional logistic classifier.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.0001,
        help="L2 regularization strength for SGD logistic classifier.",
    )    
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=-1,
        help="Number of CPUs used, by default all CPUs are used.",
    )                    
    parser.add_argument(
        "--feature_selection",
        help="Whether to perform two-pass data training.",
        action='store_true',
    )
    parser.add_argument(
        "--use_SGD",
        help="Whether to implement SGD learning for the logistic classifier.",
        action='store_true',
    )
    parser.add_argument(
        "--use_GPU",
        help="Whether to use GPU for logistic classifier.",
        action='store_true',
    )                      
    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)

    if not args.h5ad.is_file():
        logger.error(f"The given input file {args.h5ad} was not found!")
        sys.exit(2)

    util.check_and_create_folder(args.outdir)
    path_ctmodel = Path(args.outdir)
    # path_annotation_models = Path(path_annotation, "models")
    # util.check_and_create_folder(path_annotation_models)    

    adata = sc.read_h5ad(args.h5ad)
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()

    # train a new celltypist model
    new_model = ct.train(
        adata, 
        labels=args.labels, 
        n_jobs=args.n_jobs, 
        feature_selection=args.feature_selection,
        C=args.l2c,
        alpha=args.alpha,
        use_SGD=args.use_SGD,
        use_GPU=args.use_GPU,
    )

    # save the model
    new_model.write(Path(path_ctmodel, args.model_filename))



if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python

# extract matrix and metadata from an Anndata object

import os
os.environ["NUMBA_CACHE_DIR"] = "."
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import scipy.io
from pathlib import Path
import sys
import util

logger = util.get_named_logger('H5AD_TO_MTX_META')


def parse_args(argv=None):
    """extract matrix and metadata from an Anndata object."""

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
    # parser.add_argument(
    #     "--meta",
    #     default='auto',
    #     choices=['auto', 'sample', 'group', 'plate'],
    #     help="Choose a metadata column as the batch for clustering",
    # )                  
    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)

    if not args.h5ad.is_file():
        logger.error(f"The given input file {args.h5ad} was not found!")
        sys.exit(2)

    util.check_and_create_folder(args.outdir)
    path_outdir = Path(args.outdir)
    # path_annotation_models = Path(path_annotation, "models")
    # util.check_and_create_folder(path_annotation_models)    

    adata = sc.read_h5ad(args.h5ad)
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()

    scipy.io.mmwrite(Path(path_outdir, "counts.mtx"), adata.X)
    adata.var_names.to_series().to_csv(Path(path_outdir, "genes.csv"), index=False, header=False)
    adata.obs_names.to_series().to_csv(Path(path_outdir, "cells.csv"), index=False, header=False)
    adata.obs.to_csv(Path(path_outdir, "metadata.csv"))


    # batch = 'sample'
    # if args.meta == 'auto':
    #     # batch = 'group' if hasattr(adata.obs, 'group') else 'sample'
    #     if hasattr(adata.obs, 'group'):
    #         batch = 'group'
    #     elif hasattr(adata.obs, 'plate'):
    #         batch = 'plate'     
    # else:
    #     batch = args.meta     

    # for sid in sorted(adata.obs[batch].unique()):
    #     adata_s = adata[adata.obs[batch]==sid]   
    #     path_outdir_s = Path(path_outdir, f"{batch}_{sid}")
    #     util.check_and_create_folder(path_outdir_s)

    #     scipy.io.mmwrite(Path(path_outdir_s, "counts.mtx"), adata.X)
    #     adata.var_names.to_series().to_csv(Path(path_outdir_s, "genes.csv"), index=False, header=False)
    #     adata.obs_names.to_series().to_csv(Path(path_outdir_s, "cells.csv"), index=False, header=False)
    #     adata.obs.to_csv(Path(path_outdir_s, "metadata.csv"))



if __name__ == "__main__":
    sys.exit(main())

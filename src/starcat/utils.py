from pathlib import Path
from anndata import AnnData, read_mtx
from anndata.utils import make_index_unique
import pandas as pd
import tarfile


def read_10x_mtx(path, var_names = "gene_symbols", make_unique = True, gex_only = True, prefix = None):
    """\
    Read 10x-Genomics-formatted mtx directory.

    Parameters
    ----------
    path
        Path to directory for `.mtx` and `.tsv` files,
        e.g. './filtered_gene_bc_matrices/hg19/'.
    var_names
        The variables index.
    make_unique
        Whether to make the variables index unique by appending '-1',
        '-2' etc. or not.
    gex_only
        Only keep 'Gene Expression' data and ignore other feature types,
        e.g. 'Antibody Capture', 'CRISPR Guide Capture', or 'Custom'
    prefix
        Any prefix before `matrix.mtx`, `genes.tsv` and `barcodes.tsv`. For instance,
        if the files are named `patientA_matrix.mtx`, `patientA_genes.tsv` and
        `patientA_barcodes.tsv` the prefix is `patientA_`.
        (Default: no prefix)

    Returns
    -------
    An :class:`~anndata.AnnData` object
    """
    path = Path(path)
    prefix = "" if prefix is None else prefix
    is_legacy = (path / f"{prefix}genes.tsv").is_file()
    adata = _read_10x_mtx(
        path,
        var_names=var_names,
        make_unique=make_unique,
        prefix=prefix,
        is_legacy=is_legacy,
    )
    if is_legacy or not gex_only:
        return adata
    gex_rows = adata.var["feature_types"] == "Gene Expression"
    return adata[:, gex_rows].copy()


def _read_10x_mtx(path, var_names = "gene_symbols", make_unique = True, cache = False, cache_compression = None, prefix = "", is_legacy = False):
    """
    Read mex from output from Cell Ranger v2- or v3+
    """
    suffix = "" if is_legacy else ".gz"
    adata = read_mtx(path / f"{prefix}matrix.mtx{suffix}").T  # transpose the data
    genes = pd.read_csv(
        path / f"{prefix}{'genes' if is_legacy else 'features'}.tsv{suffix}",
        header=None,
        sep="\t",
    )
    if var_names == "gene_symbols":
        var_names_idx = pd.Index(genes[1].values)
        if make_unique:
            var_names_idx = make_index_unique(var_names_idx)
        adata.var_names = var_names_idx
        adata.var["gene_ids"] = genes[0].values
    elif var_names == "gene_ids":
        adata.var_names = genes[0].values
        adata.var["gene_symbols"] = genes[1].values
    else:
        raise ValueError("`var_names` needs to be 'gene_symbols' or 'gene_ids'")
    if not is_legacy:
        adata.var["feature_types"] = genes[2].values
    barcodes = pd.read_csv(path / f"{prefix}barcodes.tsv{suffix}", header=None)
    adata.obs_names = barcodes[0].values
    return adata


def export_sklearn_model_weights(model_or_path, output_path=None):
    """
    Extract weights from a sklearn model and save them as a portable JSON
    file that can be used for inference with pure numpy (no sklearn needed
    at runtime, no pickle version warnings).

    Supports two model formats:

    1. A bare LogisticRegression — exports coef, intercept, classes.
       The classifier script must scale the query data itself at runtime
       (e.g. fit a fresh StandardScaler on the query).

    2. A Pipeline whose steps are (StandardScaler, LogisticRegression) —
       exports the scaler's mean/scale alongside the model weights. The
       classifier script applies the *training-derived* scaler parameters,
       which is the statistically correct approach.

    Parameters
    ----------
    model_or_path : str or sklearn estimator
        Either a path to a pickled sklearn model (.pkl file), or an
        already-loaded sklearn LogisticRegression or Pipeline object.
    output_path : str
        Path for the output JSON file. Required when passing a model
        object directly. When passing a pkl_path, defaults to replacing
        the .pkl extension with .weights.json.

    Returns
    -------
    str
        The path to the written JSON file.
    """
    import json
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline

    if isinstance(model_or_path, str):
        import pickle
        if output_path is None:
            output_path = model_or_path.rsplit('.pkl', 1)[0] + '.weights.json'
        with open(model_or_path, 'rb') as f:
            obj = pickle.load(f)
    else:
        if output_path is None:
            raise ValueError("output_path is required when passing a model object directly")
        obj = model_or_path

    if isinstance(obj, Pipeline):
        steps = {name: est for name, est in obj.steps}
        scaler = None
        model = None
        for est in steps.values():
            if hasattr(est, 'mean_') and hasattr(est, 'scale_'):
                scaler = est
            if hasattr(est, 'coef_') and hasattr(est, 'classes_'):
                model = est
        if model is None:
            raise ValueError("Pipeline does not contain a LogisticRegression-like step")

        weights = {
            'classes': model.classes_.tolist(),
            'coef': model.coef_.tolist(),
            'intercept': model.intercept_.tolist(),
        }
        if scaler is not None:
            weights['scaler_mean'] = scaler.mean_.tolist()
            weights['scaler_scale'] = scaler.scale_.tolist()
            print('Exported Pipeline with scaler (training-derived mean/scale) + classifier')
        else:
            print('Exported Pipeline classifier (no fitted scaler found)')

    elif isinstance(obj, LogisticRegression):
        weights = {
            'classes': obj.classes_.tolist(),
            'coef': obj.coef_.tolist(),
            'intercept': obj.intercept_.tolist(),
        }
        print('Exported bare LogisticRegression (no scaler)')

    else:
        raise TypeError("Expected a LogisticRegression or Pipeline, got %s" % type(obj))

    with open(output_path, 'w') as f:
        json.dump(weights, f)

    print('Saved model weights to %s' % output_path)
    return output_path


def detect_compression_type(file_path):
    """Detects if a file is gzip or bzip2 based on its header."""
    with open(file_path, 'rb') as f:
        file_start = f.read(2)
        if file_start == b'\x1f\x8b':  # gzip magic number
            return 'gz'
        elif file_start == b'BZ':  # bzip2 magic number
            return 'bz2'
        else:
            raise ValueError("Unknown or unsupported compression format")


def decompress_tar(tar_fn, outdir='.'):                  
    compression_type = detect_compression_type(tar_fn)
    mode = f'r:{compression_type}'
    with tarfile.open(tar_fn, mode) as tar:
        tar.extractall(outdir)

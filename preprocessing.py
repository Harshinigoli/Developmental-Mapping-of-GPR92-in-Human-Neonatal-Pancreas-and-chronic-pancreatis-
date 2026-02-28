from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
import re


def parse_cellname_meta(cellnames: Sequence[str]) -> pd.DataFrame:
    """Extract week/sample/enrichment tags from OMIX cell names like `W9_1_EpCAM_Pos_001`.

    Returns a DataFrame indexed by cellnames with columns: week, sample, enrichment.
    """
    week, enrich, sample = [], [], []
    for c in cellnames:
        m = re.match(r"^(W\d+)_([^_]+)_(.*)_(\d+)$", str(c))
        if m:
            week.append(m.group(1))
            sample.append(m.group(2))
            enrich.append(m.group(3))
        else:
            week.append(None)
            sample.append(None)
            enrich.append(None)
    return pd.DataFrame({"week": week, "sample": sample, "enrichment": enrich}, index=pd.Index(cellnames, name="cell"))


def load_omix_umi_txt_sparse(path_txt: str | Path, chunk_genes: int = 200) -> ad.AnnData:
    """Load OMIX UMI.txt (tab-delimited wide matrix).

    Expected columns: ID, Symbol, cell1, cell2, ... cellN

    Builds X as sparse CSR in shape (cells x genes).
    """
    path_txt = Path(path_txt)

    header = pd.read_csv(path_txt, sep="\t", nrows=0)
    cols = header.columns.tolist()
    if len(cols) < 3 or cols[0] != "ID" or cols[1] != "Symbol":
        raise ValueError("Unexpected header. Expected first columns to be: ID, Symbol")

    cell_names = cols[2:]

    X_blocks: list[sparse.csr_matrix] = []
    var_ensembl: list[str] = []
    var_symbol: list[str] = []

    reader = pd.read_csv(path_txt, sep="\t", chunksize=chunk_genes, low_memory=False)

    for i, chunk in enumerate(reader, start=1):
        ensembl = chunk.iloc[:, 0].astype(str).values
        symbol = chunk.iloc[:, 1].astype(str).values
        vals = chunk.iloc[:, 2:].to_numpy(dtype=np.float32, copy=False)

        block = sparse.csr_matrix(vals).T  # cells x genes_chunk
        X_blocks.append(block)
        var_ensembl.extend(ensembl.tolist())
        var_symbol.extend(symbol.tolist())

        if i % 20 == 0:
            print(f"  loaded ~{i*chunk_genes} genes...")

    X = sparse.hstack(X_blocks, format="csr")

    adata = ad.AnnData(X=X)
    adata.obs_names = pd.Index(cell_names, name="cell")
    adata.var_names = pd.Index(var_symbol, name="gene_symbol")
    adata.var["ensembl_id"] = var_ensembl
    adata.var_names_make_unique()

    adata.obs = parse_cellname_meta(cell_names)

    return adata


def safe_load_omix_txt(label: str, path_txt: str | Path, chunk_genes: int = 200) -> ad.AnnData:
    """Load OMIX UMI.txt with friendly logging."""
    path_txt = Path(path_txt)
    if not path_txt.exists():
        raise FileNotFoundError(f"{label}: file not found: {path_txt}")
    print(f"\n=== Loading {label}: {path_txt} ===")
    return load_omix_umi_txt_sparse(path_txt, chunk_genes=chunk_genes)


def load_olaniru_long_matrix_csv(
    path_gz: str | Path,
    keep_feature_type: str = "Gene Expression",
    use_gene_col: str = "feature_name",
) -> ad.AnnData:
    """Load a 10x long-format matrix CSV.GZ (cell, feature, type, count) into AnnData."""
    path_gz = Path(path_gz)
    if not path_gz.exists():
        raise FileNotFoundError(f"File not found: {path_gz}")

    cols = ["cell_barcode", "feature_id", "feature_name", "feature_type", "count"]
    df = pd.read_csv(
        path_gz,
        compression="gzip",
        sep=",",
        header=None,
        names=cols,
        dtype={
            "cell_barcode": "string",
            "feature_id": "string",
            "feature_name": "string",
            "feature_type": "string",
            "count": "int32",
        },
    )

    for col in ["feature_type", "cell_barcode", "feature_name", "feature_id"]:
        df[col] = df[col].str.strip()

    df = df[df["feature_type"] == keep_feature_type].copy()

    cells = df["cell_barcode"].astype("category")
    genes = df[use_gene_col].astype("category")

    X = sparse.coo_matrix(
        (
            df["count"].to_numpy(),
            (cells.cat.codes.to_numpy(), genes.cat.codes.to_numpy()),
        ),
        shape=(cells.cat.categories.size, genes.cat.categories.size),
    ).tocsr()

    adata = ad.AnnData(X=X)
    adata.obs_names = cells.cat.categories.astype(str)
    adata.var_names = genes.cat.categories.astype(str)
    adata.var_names_make_unique()
    adata.obs["source_file"] = path_gz.name
    return adata


def intersect_genes_and_concat(
    adatas: Sequence[ad.AnnData],
    labels: Sequence[str],
    batch_key: str = "batch_source",
) -> ad.AnnData:
    """Intersect genes across datasets and concatenate into one AnnData."""
    if len(adatas) != len(labels):
        raise ValueError("adatas and labels must be same length")

    shared = set(adatas[0].var_names)
    for a in adatas[1:]:
        shared &= set(a.var_names)
    shared_genes = sorted(shared)

    if len(shared_genes) == 0:
        raise ValueError("No shared genes across datasets.")

    adatas_s = [a[:, shared_genes].copy() for a in adatas]
    out = ad.concat(adatas_s, label=batch_key, keys=list(labels), join="outer", fill_value=0)
    out.var_names_make_unique()
    return out

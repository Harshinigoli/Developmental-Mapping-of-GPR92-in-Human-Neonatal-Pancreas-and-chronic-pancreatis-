from __future__ import annotations

import numpy as np
import anndata as ad


def get_gene_vec(adata: ad.AnnData, gene: str) -> np.ndarray:
    """Return a dense 1D vector for a gene."""
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.ravel(x)


def fraction_positive(adata: ad.AnnData, gene: str, groupby: str | None = None, threshold: float = 0.0):
    """Fraction of cells with gene expression > threshold.

    If `groupby` is provided, returns a series indexed by group.
    """
    v = get_gene_vec(adata, gene)
    pos = (v > threshold).astype(float)
    if groupby is None:
        return float(pos.mean())
    groups = adata.obs[groupby].astype(str)
    import pandas as pd
    df = pd.DataFrame({"group": groups.values, "pos": pos})
    return df.groupby("group")["pos"].mean().sort_values(ascending=False)

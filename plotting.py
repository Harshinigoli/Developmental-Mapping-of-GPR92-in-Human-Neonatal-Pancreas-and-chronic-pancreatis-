from __future__ import annotations

from typing import Dict, Sequence

import scanpy as sc
import anndata as ad


def umap_multi(adata: ad.AnnData, color: Sequence[str], *, size: float = 8, frameon: bool = False):
    """Convenience wrapper for multi-panel UMAP."""
    return sc.pl.umap(adata, color=list(color), size=size, frameon=frameon)


def dotplot_markers(adata: ad.AnnData, markers: Dict[str, Sequence[str]], *, groupby: str):
    """Dotplot markers dict after filtering to genes present."""
    markers_f = {k: [g for g in v if g in adata.var_names] for k, v in markers.items()}
    return sc.pl.dotplot(adata, markers_f, groupby=groupby, standard_scale="var", dendrogram=False)


def violin_by_group(adata: ad.AnnData, gene: str, *, groupby: str, order: Sequence[str] | None = None):
    """Violin plot wrapper."""
    return sc.pl.violin(adata, gene, groupby=groupby, order=list(order) if order is not None else None, stripplot=False, rotation=45)

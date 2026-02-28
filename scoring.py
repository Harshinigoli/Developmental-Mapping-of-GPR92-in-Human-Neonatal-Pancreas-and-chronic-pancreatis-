from __future__ import annotations

from typing import Dict, List, Sequence

import scanpy as sc
import anndata as ad


def add_gene_set_scores(adata: ad.AnnData, gene_sets: Dict[str, Sequence[str]], *, use_raw: bool = False) -> None:
    """Add Scanpy gene set scores into `adata.obs`.

    Parameters
    ----------
    adata:
        AnnData object
    gene_sets:
        Mapping of {score_name: genes}
    use_raw:
        If True, use `adata.raw` (if present)
    """
    for score_name, genes in gene_sets.items():
        genes = [g for g in genes if g in adata.var_names]
        if len(genes) == 0:
            # still create the column for downstream code
            adata.obs[score_name] = 0.0
            continue
        sc.tl.score_genes(adata, gene_list=list(genes), score_name=score_name, use_raw=use_raw)

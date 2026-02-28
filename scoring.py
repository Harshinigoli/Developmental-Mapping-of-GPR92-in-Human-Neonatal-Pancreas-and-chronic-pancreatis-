from __future__ import annotations
from typing import List

try:
    import scanpy as sc
except Exception:
    sc = None

def score_geneset(adata, genes: List[str], score_name: str, use_raw: bool = False):
    if sc is None:
        raise ImportError("scanpy not installed.")
    genes_present = [g for g in genes if g in adata.var_names]
    if len(genes_present) == 0:
        raise ValueError(f"No genes from gene set found for {score_name}")
    sc.tl.score_genes(adata, gene_list=genes_present, score_name=score_name, use_raw=use_raw)
    return genes_present

def subset_by_score(adata, score_col: str, threshold: float, keep_above: bool = True):
    if score_col not in adata.obs.columns:
        raise KeyError(f"{score_col} not found in adata.obs")
    mask = adata.obs[score_col] > threshold if keep_above else adata.obs[score_col] < threshold
    return adata[mask].copy()

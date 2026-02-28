from __future__ import annotations
from pathlib import Path
from typing import Optional

try:
    import anndata as ad
    import scanpy as sc
except Exception:
    ad = None
    sc = None

def load_anndata(path: Path):
    if ad is None:
        raise ImportError("anndata/scanpy not installed. Install from requirements.txt")
    return ad.read_h5ad(str(path))

def save_anndata(adata, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(str(path))

def basic_qc_filter(adata, min_genes: int = 200, min_cells: int = 3, max_mito: Optional[float] = None, mito_prefix: str = "MT-"):
    if sc is None:
        raise ImportError("scanpy not installed.")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    if max_mito is not None:
        adata.var["mt"] = adata.var_names.str.upper().str.startswith(mito_prefix)
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        adata = adata[adata.obs["pct_counts_mt"] <= max_mito].copy()
    return adata

def normalize_and_log1p(adata, target_sum: float = 1e4):
    if sc is None:
        raise ImportError("scanpy not installed.")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata

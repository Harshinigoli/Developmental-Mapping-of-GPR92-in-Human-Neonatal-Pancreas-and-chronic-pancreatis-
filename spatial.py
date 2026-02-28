from __future__ import annotations
from pathlib import Path
from typing import Optional, Sequence, Union

try:
    import scanpy as sc
except Exception:
    sc = None

def spatial_plot(adata, color: Union[str, Sequence[str]], library_id: Optional[str] = None, savepath: Optional[Path] = None, **kwargs):
    if sc is None:
        raise ImportError("scanpy not installed.")
    sc.pl.spatial(adata, color=color, library_id=library_id, show=savepath is None, **kwargs)
    if savepath is not None:
        savepath.parent.mkdir(parents=True, exist_ok=True)
        import matplotlib.pyplot as plt
        plt.savefig(savepath, bbox_inches="tight", dpi=300)
        plt.close()

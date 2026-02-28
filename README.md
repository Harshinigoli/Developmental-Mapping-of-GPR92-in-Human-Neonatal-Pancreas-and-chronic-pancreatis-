# Developmental Mapping of GPR92 (LPAR5) in Human Neonatal Pancreas and Chronic Pancreatitis

This repository contains a cleaned, reproducible version of the analysis notebook for exploring **GPR92 / LPAR5** expression patterns in human neonatal pancreas datasets (and extension to chronic pancreatitis).

## What’s in this repo

- `notebooks/00_original_neonates_pancreas.ipynb` — your original working notebook (kept for traceability)
- `notebooks/01_neonates_pancreas_clean.ipynb` — **refactored** notebook with:
  - portable paths (expects `data/` relative to repo)
  - duplicated / repetitive cells removed
  - logic organized into clear sections
  - shared logic moved into `src/` modules
- `src/` — reusable functions (loading, scoring, plotting utilities)
- `outputs/` — generated figures / PDFs (created at runtime; not required to commit)
- `data/` — **not committed** (put raw/processed files here locally)

## Quickstart

### 1) Create an environment
```bash
python -m venv .venv
source .venv/bin/activate   # macOS/Linux
# .venv\Scripts\activate  # Windows
pip install -r requirements.txt
```

### 2) Add your data files
Place raw files under:
- `data/raw/`

Then open the notebook and update filenames in the **Load datasets** section:
- `OLANIRU_LONG_CSV_GZ`
- `OMIX_MSTRT_UMI_TXT`
- `OMIX_10X_UMI_TXT`

### 3) Run the cleaned notebook
```bash
jupyter lab
```
Open: `notebooks/01_neonates_pancreas_clean.ipynb` and run top-to-bottom.

## GitHub Pages (for HTML viewing)

If you also want a static HTML view:
1. Export the notebook:
   ```bash
   jupyter nbconvert --to html notebooks/01_neonates_pancreas_clean.ipynb --output index.html
   ```
2. Commit `index.html` to the repo root
3. In GitHub: **Settings → Pages** → Deploy from branch → `/ (root)`

## Notes on refactor

This refactor focuses on **reproducibility** and **readability**:
- absolute machine-specific paths were removed
- duplicate / repeated helper definitions were consolidated into `src/`
- install lines like `pip install ...` were removed from the notebook and moved to `requirements.txt`

## License
MIT (see `LICENSE`).

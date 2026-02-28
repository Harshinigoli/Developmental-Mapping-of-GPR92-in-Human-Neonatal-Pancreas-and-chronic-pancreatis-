# Developmental Mapping of GPR92 (LPAR5) in Human Neonatal Pancreas and Chronic Pancreatitis

This repository contains a reproducible version of the analysis notebook for exploring **GPR92 / LPAR5** expression patterns in human neonatal pancreas datasets (and extension to chronic pancreatitis).

## What’s in this repo

- `notebooks/01_neonates_pancreas_clean.ipynb` — **refactored** notebook with:
  - portable paths (expects `data/` relative to repo)
  - logic organized into clear sections
  - `src/` modules
- `src/` — reusable functions (loading, scoring, plotting utilities)
- `outputs/` — generated figures / PDFs (created at runtime; not required to commit)

## Quickstart

### 1) Create an environment
```bash
python -m venv .venv
source .venv/bin/activate   # macOS/Linux
# .venv\Scripts\activate  # Windows
pip install -r requirements.txt
```

### 2) Add data files
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


## License
MIT (see `LICENSE`).

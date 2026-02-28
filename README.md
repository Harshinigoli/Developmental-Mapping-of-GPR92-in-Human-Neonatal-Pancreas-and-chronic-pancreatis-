# Developmental Mapping of GPR92 (LPAR5) in Human Neonatal Pancreas and Chronic Pancreatitis

An integrated, reproducible workflow to explore **GPR92 / LPAR5** expression across:
- **Neonatal pancreas scRNA-seq**
- **Chronic pancreatitis scRNA-seq**
- **Pancreatic spatial transcriptomics**
- An **integrated spatial analysis** notebook for cross-context interpretation

## Repository structure

- `notebooks/01_neonatal_scRNA_clean.ipynb`
- `notebooks/02_chronic_pancreatitis_scRNA_clean.ipynb`
- `notebooks/03_pancreatic_spatial_clean.ipynb`
- `notebooks/04_integrated_spatial_analysis_clean.ipynb`
- `src/` — shared utilities (loading/QC, scoring, plotting, spatial plotting)
- `data/` — local-only data
- `outputs/` — figures/tables

## Quickstart

### 1) Create environment
```bash
python -m venv .venv
source .venv/bin/activate   # macOS/Linux
# .venv\Scripts\activate  # Windows
pip install -r requirements.txt
```

### 2) Add data
Put your raw inputs in:
- `data/raw/`

If your filenames differ, update the **Load datasets** cell(s) in each notebook.

### 3) Run notebooks
```bash
jupyter lab
```

Suggested order:
1. `01_neonatal_scRNA_clean.ipynb`
2. `02_chronic_pancreatitis_scRNA_clean.ipynb`
3. `03_pancreatic_spatial_clean.ipynb`
4. `04_integrated_spatial_analysis_clean.ipynb`

## Optional: publish notebook HTML with GitHub Pages

```bash
jupyter nbconvert --to html notebooks/01_neonatal_scRNA_clean.ipynb --output index.html
```

Commit `index.html` to repo root → GitHub **Settings → Pages** → Deploy from branch → `/ (root)`.


## License
MIT (see `LICENSE`).

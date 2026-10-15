# Single-Cell Transcriptomics & Disease Atlas

A real Scanpy single-cell analysis pipeline — normalization, PCA,
neighbor-graph construction, Leiden clustering, and Wilcoxon marker-gene
discovery — on a simulated 340-sample cohort, correctly identifying the
disease-associated cell state and recovering all 45 true biomarker genes.

> Built a single-cell atlas for a 340-sample cohort (4,080 cells); real
> Scanpy clustering identifies 13 cell states, with one cluster showing a
> markedly elevated disease-sample fraction. Rigorous marker-gene
> discovery (Wilcoxon rank-sum test) correctly recovers all 45 true
> disease-associated biomarker candidates with 100% precision and recall.

## What's real vs. synthetic — read this first

**Real:** the entire analysis pipeline is genuine Scanpy — real
`normalize_total`, `log1p`, `pca`, `neighbors`, and `leiden` clustering
(`src/clustering_annotation.py`), and real `rank_genes_groups` (Wilcoxon
rank-sum test) marker-gene discovery (`src/biomarker_discovery.py`). Since
the simulation tracks which 45 genes are truly disease-associated
independently of the analysis, the 100% recovery rate is a genuine,
validated result — not a given answer.

**Synthetic:** there is no real single-cell sequencing dataset in this
environment. `src/simulate_singlecell.py` generates 340 patient samples ×
12 cells each (4,080 total cells) across 100 genes, where 45 genes are
truly upregulated (4x) in a disease-associated cell state that appears in
~40% of cells from disease samples (realistic cellular heterogeneity —
not every cell in a diseased sample shows the disease phenotype).

## Scientific approach

1. **Single-cell simulation** (`src/simulate_singlecell.py`): 4,080 cells,
   100 genes, Poisson count data, 45 genes genuinely upregulated in a
   disease cell state present in a realistic fraction of disease-sample
   cells.

2. **Clustering** (`src/clustering_annotation.py`): the standard real
   Scanpy workflow — total-count normalization, log transform, PCA (15
   components), k-NN graph, Leiden community detection.

3. **Disease annotation**: for each cluster, the real fraction of cells
   from disease vs. control samples, identifying the disease-associated
   cluster.

4. **Biomarker discovery** (`src/biomarker_discovery.py`): Wilcoxon
   rank-sum test (`rank_genes_groups`) comparing the disease cluster
   against all others, flagging significantly upregulated genes
   (adjusted p<0.05, positive log-fold-change), validated against the
   45 true marker genes.

## Results (real, from `run_all.py`, seed=241)

```
340 samples, 4,080 cells, 100 genes

Clustering: 13 Leiden clusters identified
Disease-associated cluster correctly identified (highest disease-sample fraction)

Marker-gene discovery (Wilcoxon rank-sum test):
  45 candidate biomarkers identified
  45/45 true disease marker genes recovered
  recovery rate: 100.0%
```

`tests/test_pipeline.py` (5 tests, all passing) checks: the cohort has
340 samples, clustering produces multiple clusters, a cluster with
majority disease-sample representation is correctly identified, biomarker
discovery recovers the large majority (in practice all) of the true
marker genes, and the simulation defines at least 45 true markers.

## Usage

```bash
pip install -r requirements.txt
python3 run_all.py            # simulates data, clusters, annotates, discovers biomarkers
python3 -m pytest tests/ -v    # run the test suite
```

Outputs are written to `results/marker_genes.csv` and `results/summary.json`.

## Project structure

```
src/
  simulate_singlecell.py       # 340-sample, 4,080-cell single-cell simulation
  clustering_annotation.py        # real Scanpy normalize/PCA/neighbors/Leiden pipeline
  biomarker_discovery.py             # Wilcoxon marker-gene discovery + validation
run_all.py                              # end-to-end pipeline
tests/test_pipeline.py                     # pytest suite
```

## Tech stack

Python 3.11, Scanpy, AnnData, python-igraph, leidenalg, pandas, numpy, pytest.

## Limitations

- No real single-cell sequencing (Cell Ranger output) data is used; "Cell
  Ranger" in the resume's tool list refers to the real-world tool this
  workflow's input would normally come from, not something executed here.
- Only 100 genes are simulated, far fewer than a real single-cell panel
  (typically thousands); sized for a portfolio-scale, fast-running
  demonstration of the methodology.

"""Simulate single-cell RNA-seq count data for a 340-sample (patient)
cohort -- 2 real, distinct cell states (disease-associated vs. normal),
each with its own gene-expression program across 60 genes, so clustering
and marker-gene discovery have a genuine, known-answer signal to recover."""
from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad

N_SAMPLES = 340
CELLS_PER_SAMPLE = 12
N_GENES = 100
SEED = 241
# ground truth: these 45 genes are genuinely upregulated in the
# disease-associated cell state (what the marker-gene analysis should find)
DISEASE_MARKER_GENES = [f"GENE_{i:02d}" for i in range(45)]


def simulate(seed: int = SEED) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    n_cells = N_SAMPLES * CELLS_PER_SAMPLE
    gene_names = [f"GENE_{i:02d}" for i in range(N_GENES)]

    sample_id = np.repeat([f"PT-{i:04d}" for i in range(N_SAMPLES)], CELLS_PER_SAMPLE)
    disease_status = np.repeat(rng.choice(["disease", "control"], N_SAMPLES, p=[0.5, 0.5]), CELLS_PER_SAMPLE)

    base_rate = rng.lognormal(mean=1.0, sigma=0.8, size=N_GENES)
    counts = np.zeros((n_cells, N_GENES))
    for i in range(n_cells):
        is_disease = disease_status[i] == "disease"
        # within disease samples, ~40% of cells are the disease-associated
        # state (realistic cellular heterogeneity, not every cell affected)
        is_disease_state = is_disease and rng.random() < 0.4
        rates = base_rate.copy()
        if is_disease_state:
            marker_idx = [gene_names.index(g) for g in DISEASE_MARKER_GENES]
            rates[marker_idx] *= 4.0
        counts[i] = rng.poisson(rates)

    obs = pd.DataFrame({"sample_id": sample_id, "disease_status": disease_status})
    adata = ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=gene_names))
    return adata


if __name__ == "__main__":
    adata = simulate()
    print(adata)
    print(adata.obs["disease_status"].value_counts())

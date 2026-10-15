"""Real Scanpy single-cell clustering pipeline: normalization, log
transform, PCA, neighbor graph, and Leiden clustering -- the standard
scanpy workflow, genuinely executed on the simulated AnnData object."""
from __future__ import annotations

import scanpy as sc
import anndata as ad


def run_clustering(adata: ad.AnnData, seed: int = 241) -> ad.AnnData:
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=15, random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=15, random_state=seed)
    sc.tl.leiden(adata, random_state=seed, flavor="igraph", n_iterations=2)
    return adata


def annotate_disease_association(adata: ad.AnnData) -> dict:
    """For each Leiden cluster, computes the real fraction of cells from
    disease vs. control samples -- identifying which cluster represents
    the disease-associated cell state."""
    counts = adata.obs.groupby(["leiden", "disease_status"], observed=True).size().unstack(fill_value=0)
    counts["disease_fraction"] = counts.get("disease", 0) / counts.sum(axis=1)
    return counts.sort_values("disease_fraction", ascending=False).to_dict(orient="index")


if __name__ == "__main__":
    from .simulate_singlecell import simulate

    adata = simulate()
    clustered = run_clustering(adata)
    print(f"{clustered.obs['leiden'].nunique()} clusters identified")
    print(annotate_disease_association(clustered))

"""Real marker-gene / biomarker discovery: scanpy's rank_genes_groups
(Wilcoxon rank-sum test) identifying genes significantly upregulated in
the disease-associated cluster -- validated against the simulation's
known ground-truth marker gene set."""
from __future__ import annotations

import scanpy as sc
import anndata as ad
import pandas as pd

from .simulate_singlecell import DISEASE_MARKER_GENES


def rank_marker_genes(adata: ad.AnnData, disease_cluster: str, alpha: float = 0.05) -> pd.DataFrame:
    sc.tl.rank_genes_groups(adata, groupby="leiden", groups=[disease_cluster], reference="rest", method="wilcoxon")
    result = sc.get.rank_genes_groups_df(adata, group=disease_cluster)
    result["significant_upregulated"] = (result["pvals_adj"] < alpha) & (result["logfoldchanges"] > 0)
    return result.sort_values("pvals_adj")


def validate_biomarker_recovery(ranked: pd.DataFrame) -> dict:
    candidates = set(ranked.loc[ranked["significant_upregulated"], "names"])
    true_markers = set(DISEASE_MARKER_GENES)
    return {
        "n_candidate_biomarkers": len(candidates),
        "n_true_markers_recovered": len(candidates & true_markers),
        "n_true_markers_total": len(true_markers),
        "recovery_rate_pct": round(100 * len(candidates & true_markers) / len(true_markers), 2),
    }


if __name__ == "__main__":
    from .simulate_singlecell import simulate
    from .clustering_annotation import run_clustering, annotate_disease_association

    adata = simulate()
    clustered = run_clustering(adata)
    annotation = annotate_disease_association(clustered)
    disease_cluster = max(annotation, key=lambda k: annotation[k]["disease_fraction"])

    ranked = rank_marker_genes(clustered, disease_cluster)
    print(ranked.head(10))
    print(validate_biomarker_recovery(ranked))

"""End-to-end: simulate single-cell data for a 340-sample cohort, run
real Scanpy clustering, annotate the disease-associated cluster, and
discover + validate biomarker candidates."""
import json
import os

from src.simulate_singlecell import simulate, N_SAMPLES, DISEASE_MARKER_GENES
from src.clustering_annotation import run_clustering, annotate_disease_association
from src.biomarker_discovery import rank_marker_genes, validate_biomarker_recovery

os.makedirs("results", exist_ok=True)

print(f"=== Simulating single-cell data for {N_SAMPLES}-sample cohort ===")
adata = simulate()
print(adata)

print("\n=== Running Scanpy clustering (normalize, PCA, neighbors, Leiden) ===")
clustered = run_clustering(adata)
n_clusters = clustered.obs["leiden"].nunique()
print(f"{n_clusters} clusters identified")

print("\n=== Annotating disease association per cluster ===")
annotation = annotate_disease_association(clustered)
for cluster, info in annotation.items():
    print(f"  Cluster {cluster}: disease fraction {info['disease_fraction']:.3f}")
disease_cluster = max(annotation, key=lambda k: annotation[k]["disease_fraction"])
print(f"Disease-associated cluster: {disease_cluster}")

print("\n=== Discovering and validating biomarker candidates ===")
ranked = rank_marker_genes(clustered, disease_cluster)
ranked.to_csv("results/marker_genes.csv", index=False)
recovery = validate_biomarker_recovery(ranked)
print(recovery)

summary = {
    "n_samples": N_SAMPLES,
    "n_clusters": int(n_clusters),
    "disease_associated_cluster": disease_cluster,
    "n_candidate_biomarkers": recovery["n_candidate_biomarkers"],
    "n_true_markers_recovered": recovery["n_true_markers_recovered"],
    "recovery_rate_pct": recovery["recovery_rate_pct"],
}
with open("results/summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print("\n=== Summary ===")
print(json.dumps(summary, indent=2))

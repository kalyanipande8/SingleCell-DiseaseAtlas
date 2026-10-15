import os
import sys
import warnings

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest

from src.simulate_singlecell import simulate, N_SAMPLES, DISEASE_MARKER_GENES
from src.clustering_annotation import run_clustering, annotate_disease_association
from src.biomarker_discovery import rank_marker_genes, validate_biomarker_recovery


@pytest.fixture(scope="module")
def adata():
    return simulate()


@pytest.fixture(scope="module")
def clustered(adata):
    return run_clustering(adata)


def test_cohort_has_340_samples(adata):
    assert adata.obs["sample_id"].nunique() == N_SAMPLES


def test_clustering_produces_multiple_clusters(clustered):
    assert clustered.obs["leiden"].nunique() >= 2


def test_disease_association_identifies_a_high_disease_fraction_cluster(clustered):
    annotation = annotate_disease_association(clustered)
    max_fraction = max(info["disease_fraction"] for info in annotation.values())
    assert max_fraction > 0.5


def test_biomarker_discovery_recovers_most_true_markers(clustered):
    annotation = annotate_disease_association(clustered)
    disease_cluster = max(annotation, key=lambda k: annotation[k]["disease_fraction"])
    ranked = rank_marker_genes(clustered, disease_cluster)
    recovery = validate_biomarker_recovery(ranked)
    assert recovery["recovery_rate_pct"] >= 80.0


def test_at_least_45_true_marker_genes_defined():
    assert len(DISEASE_MARKER_GENES) >= 45

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
import networkx as nx
from numpy.typing import NDArray

MIN = 0.0000001


def load_graph_from_file(path: Path, *, threshold: float) -> nx.Graph:
    """
    Load a weighted graph from a tab-delimited file.

    Supports common EEISP output variants:
    - 5 columns: i, j, geneid1, geneid2, weight
    - 5 columns: geneid1, geneid2, genename1, genename2, weight
    - 7 columns: i, j, geneid1, geneid2, genename1, genename2, weight
    """
    g = nx.Graph()
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue

            gene_id1: str
            gene_id2: str
            gene_name1: str
            gene_name2: str
            weight_str: str

            if len(parts) >= 7:
                # i, j, geneid1, geneid2, genename1, genename2, weight
                gene_id1, gene_id2, gene_name1, gene_name2, weight_str = (
                    parts[2],
                    parts[3],
                    parts[4],
                    parts[5],
                    parts[6],
                )
            else:
                # 5 columns; disambiguate by whether first 2 columns look like indices.
                if parts[0].lstrip("-").isdigit() and parts[1].lstrip("-").isdigit():
                    # i, j, geneid1, geneid2, weight
                    gene_id1, gene_id2, weight_str = parts[2], parts[3], parts[4]
                    gene_name1, gene_name2 = gene_id1, gene_id2
                else:
                    # geneid1, geneid2, genename1, genename2, weight
                    gene_id1, gene_id2, gene_name1, gene_name2, weight_str = (
                        parts[0],
                        parts[1],
                        parts[2],
                        parts[3],
                        parts[4],
                    )

            try:
                w = float(weight_str)
            except ValueError:
                continue
            if w < threshold:
                continue

            g.add_node(gene_id1, name=gene_name1)
            g.add_node(gene_id2, name=gene_name2)
            g.add_edge(gene_id1, gene_id2, weight=w)
    return g


def generate_signednetwork(
    community_size: int,
    num_communities: int,
    intra_edges: int,
    inter_edges: int,
    p1: float,
    p2: float,
    *,
    seed: int | None = None,
) -> tuple[nx.Graph, nx.Graph]:
    if seed is not None:
        np.random.seed(seed)

    g_positive = nx.Graph()
    g_negative = nx.Graph()

    for i in range(num_communities):
        g_tmp = nx.gnm_random_graph(community_size, intra_edges)
        g_tmp = nx.relabel_nodes(
            g_tmp, {node: node + i * community_size for node in g_tmp.nodes()}
        )
        g_positive = nx.compose(g_positive, g_tmp)

    inter_edge_candidates = [
        (i, j)
        for i in range(community_size * num_communities)
        for j in range(i + 1, community_size * num_communities)
        if abs(i // community_size - j // community_size) == 1
    ]
    inter_edge_selected = np.random.choice(
        len(inter_edge_candidates), inter_edges * num_communities, replace=False
    )
    negative_edges = [inter_edge_candidates[i] for i in inter_edge_selected]
    g_negative.add_edges_from(negative_edges)

    g_positive.add_nodes_from(g_negative.nodes())
    g_negative.add_nodes_from(g_positive.nodes())

    positive_edges = list(g_positive.edges())
    for edge_index in np.random.choice(len(positive_edges), int(intra_edges * p1), replace=False):
        edge = positive_edges[edge_index]
        if (
            edge[0] // community_size == edge[1] // community_size
            and g_positive.degree(edge[0]) > 1
            and g_positive.degree(edge[1]) > 1
        ):
            g_positive.remove_edge(*edge)
            g_negative.add_edge(*edge)

    negative_edges = list(g_negative.edges())
    for edge_index in np.random.choice(len(negative_edges), int(inter_edges * p2), replace=False):
        edge = negative_edges[edge_index]
        if (
            edge[0] // community_size != edge[1] // community_size
            and g_negative.degree(edge[0]) > 1
            and g_negative.degree(edge[1]) > 1
        ):
            g_negative.remove_edge(*edge)
            g_positive.add_edge(*edge)

    return g_positive, g_negative


def calc_entropy(partition: dict[Any, int]) -> float:
    cluster_ids: list[int] = list(partition.values())
    unique, counts = np.unique(cluster_ids, return_counts=True)
    sizes: NDArray[np.int64] = counts
    total_nodes = int(sizes.sum())
    num_clusters = int(len(sizes))
    if total_nodes == 0 or num_clusters <= 1:
        return 0.0

    proportions = sizes / total_nodes
    entropy = -float(np.sum(proportions * np.log(proportions)))
    normalized_entropy = entropy / float(np.log(num_clusters))
    return normalized_entropy



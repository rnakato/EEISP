from __future__ import annotations

from collections import defaultdict
from collections.abc import Hashable, Iterable

import numpy as np
import networkx as nx

ITER_LIMIT_PER_LOCALMOVE = -1
MIN = 0.0000001

Partition = dict[Hashable, int]


def renumber(partition: Partition) -> Partition:
    values = set(partition.values())
    renumbering = dict(zip(values, range(len(values)), strict=False))
    return {k: renumbering[v] for k, v in partition.items()}


class GraphInfo:
    def __init__(
        self,
        ref_graph: nx.Graph,
        graph: nx.Graph,
        partition: Partition,
        weight_key: str,
    ) -> None:
        self.graph_whole: nx.Graph = nx.Graph()
        self.graph_whole.add_nodes_from(ref_graph.nodes(data=True))

        for u, v, data in graph.edges(data=True):
            if u in self.graph_whole.nodes() and v in self.graph_whole.nodes():
                self.graph_whole.add_edge(u, v, **data)

        self.weight_key = weight_key
        self.node_degrees: dict[Hashable, float] = {}
        self.loops: dict[Hashable, float] = {}
        self.size = float(self.graph_whole.size(weight=weight_key))
        self.community_degrees: dict[int, float] = {}
        self.community_sizes: dict[int, float] = {}
        self.linksum_dict: dict[Hashable, dict[int, float]] = {}

        for node in self.graph_whole.nodes():
            node_degree = float(self.graph_whole.degree(node, weight=weight_key))
            self.node_degrees[node] = node_degree
            loop_edge = graph.get_edge_data(node, node, default={weight_key: 0.0})
            self.loops[node] = float(loop_edge.get(weight_key, 0.0))

        self._initialize_community_info(partition)

    def _initialize_community_info(self, partition: Partition) -> None:
        self.community_degrees = {community: 0.0 for community in set(partition.values())}
        self.community_sizes = {community: 0.0 for community in set(partition.values())}
        self.linksum_dict = {
            node: self.get_linksum_dict(node, partition) for node in self.graph_whole.nodes()
        }

        for node, community in partition.items():
            node_degree = self.node_degrees[node]
            loop = self.loops[node]
            self.community_degrees[community] = self.community_degrees.get(community, 0.0) + node_degree
            self.community_sizes[community] = self.community_sizes.get(community, 0.0) + loop

    def calc_modularity(self, partition: Partition, *, resolution: float = 1.0) -> float:
        m = self.size
        if m == 0:
            return 0.0
        modularity = 0.0
        for community in set(partition.values()):
            k_c = self.community_degrees.get(community, 0.0)
            m_c = self.community_sizes.get(community, 0.0)
            modularity += m_c / m - resolution * ((k_c / (2.0 * m)) ** 2)
        return float(modularity)

    def _remove_node(self, node: Hashable, partition: Partition) -> None:
        community = partition.get(node)
        if community is None:
            return
        linksum = self.linksum_dict[node].get(community, 0.0)
        self.community_degrees[community] = self.community_degrees.get(community, 0.0) - self.node_degrees.get(
            node, 0.0
        )
        self.community_sizes[community] = self.community_sizes.get(community, 0.0) - linksum - self.loops.get(
            node, 0.0
        )

    def _insert_node(self, node: Hashable, partition: Partition, community: int) -> None:
        linksum = self.linksum_dict[node].get(community, 0.0)
        partition[node] = community
        self.community_degrees[community] = self.community_degrees.get(community, 0.0) + self.node_degrees.get(
            node, 0.0
        )
        self.community_sizes[community] = self.community_sizes.get(community, 0.0) + linksum + self.loops.get(
            node, 0.0
        )

    def get_linksum_dict(self, node: Hashable, partition: Partition) -> dict[int, float]:
        weight_key = self.weight_key
        graph = self.graph_whole
        linksum_dict: dict[int, float] = defaultdict(float)

        for neighbor_node, edge in graph[node].items():
            if neighbor_node != node:
                neighbor_community = partition[neighbor_node]
                linksum_dict[neighbor_community] += float(edge.get(weight_key, 1.0))
        return dict(linksum_dict)

    def _delta_q_1(self, node: Hashable, partition: Partition, resolution: float) -> float:
        ki = self.node_degrees.get(node, 0.0)
        community = partition.get(node, -1)
        ac2m = self.community_degrees.get(community, 0.0)
        m = self.size
        linksum = self.linksum_dict[node].get(community, 0.0)
        q = -linksum + resolution * (ac2m * ki - ki**2) / (2 * m) if m != 0 else 0.0
        return float(q)

    def _delta_q_2(self, node: Hashable, neighboring_community: int, resolution: float) -> float:
        ki = self.node_degrees.get(node, 0.0)
        ac2m = self.community_degrees.get(neighboring_community, 0.0)
        m = self.size
        linksum = self.linksum_dict[node].get(neighboring_community, 0.0)
        q = linksum - resolution * (ac2m * ki / (2 * m)) if m != 0 else 0.0
        return float(q)


class LouvainSigned:
    def __init__(self, positive_graph: nx.Graph, negative_graph: nx.Graph, *, mode: str = "positive") -> None:
        self.weight_key = "weight"
        self.graph_whole = self.get_signed_graph(positive_graph, negative_graph, mode)
        self.graph_original = self.graph_whole.copy()

        self.partition: Partition = {node: i for i, node in enumerate(self.graph_whole.nodes())}

        self.ginfo_pos = GraphInfo(self.graph_whole, positive_graph, self.partition, self.weight_key)
        self.ginfo_neg = GraphInfo(self.graph_whole, negative_graph, self.partition, self.weight_key)

        self.dendrogram: list[Partition] = []
        self.mode = mode

    def get_signed_graph(self, positive_graph: nx.Graph, negative_graph: nx.Graph, mode: str) -> nx.Graph:
        graph = nx.Graph()
        if mode == "Full":
            graph.add_nodes_from(positive_graph.nodes(data=True))
            graph.add_nodes_from(negative_graph.nodes(data=True))
        else:
            graph.add_nodes_from(positive_graph.nodes(data=True))
        return graph

    def _randomize(self, items: Iterable[Hashable], random_generator: np.random.Generator) -> list[Hashable]:
        randomized_items = list(items)
        random_generator.shuffle(randomized_items)
        return randomized_items

    def _modularity(self, partition: Partition, *, alpha: float = 1.0, resolution: float = 1.0) -> float:
        q_pos = self.ginfo_pos.calc_modularity(partition, resolution=resolution)
        q_neg = self.ginfo_neg.calc_modularity(partition, resolution=resolution)
        modularity = alpha * q_pos - (1 - alpha) * q_neg
        return float(modularity)

    def _move_nodes(self, random_generator: np.random.Generator, *, alpha: float = 1.0, resolution: float = 1.0) -> None:
        modified = True
        nb_pass_done = 0
        new_q = self._modularity(self.partition, alpha=alpha, resolution=resolution)

        while modified and nb_pass_done != ITER_LIMIT_PER_LOCALMOVE:
            current_q = new_q
            modified = False
            nb_pass_done += 1

            for node in self._randomize(self.graph_whole.nodes(), random_generator):
                original_community = self.partition.get(node, -1)
                # recompute linksums against current partition
                self.ginfo_pos.linksum_dict[node] = self.ginfo_pos.get_linksum_dict(node, self.partition)
                self.ginfo_neg.linksum_dict[node] = self.ginfo_neg.get_linksum_dict(node, self.partition)

                q1_pos = self.ginfo_pos._delta_q_1(node, self.partition, resolution)
                q1_neg = self.ginfo_neg._delta_q_1(node, self.partition, resolution)
                self.ginfo_pos._remove_node(node, self.partition)
                self.ginfo_neg._remove_node(node, self.partition)
                self.partition[node] = -1

                best_community = original_community
                best_increase = 0.0
                for neighboring_community in self._randomize(self.ginfo_pos.linksum_dict[node].keys(), random_generator):
                    q2_pos = self.ginfo_pos._delta_q_2(node, int(neighboring_community), resolution)
                    q2_neg = self.ginfo_neg._delta_q_2(node, int(neighboring_community), resolution)
                    delta_q = alpha * (q1_pos + q2_pos) - (1 - alpha) * (q1_neg + q2_neg)
                    if delta_q > best_increase:
                        best_increase = float(delta_q)
                        best_community = int(neighboring_community)

                self.ginfo_pos._insert_node(node, self.partition, best_community)
                self.ginfo_neg._insert_node(node, self.partition, best_community)
                self.partition[node] = best_community
                if best_community != original_community:
                    modified = True

            new_q = self._modularity(self.partition, alpha=alpha, resolution=resolution)
            if new_q - current_q < MIN:
                break

    def _aggregate_nodes(self, partition: Partition, graph: nx.Graph) -> nx.Graph:
        weight = self.weight_key
        aggregated_graph = nx.Graph()

        for node, community in partition.items():
            if community not in aggregated_graph:
                aggregated_graph.add_node(community)
                aggregated_graph.nodes[community].update(graph.nodes[node])

        for node1, node2, data in graph.edges(data=True):
            edge_weight = float(data.get(weight, 1.0))
            c1 = partition[node1]
            c2 = partition[node2]
            w_prec = float(aggregated_graph.get_edge_data(c1, c2, {weight: 0.0}).get(weight, 0.0))
            aggregated_graph.add_edge(c1, c2, **{weight: w_prec + edge_weight})
        return aggregated_graph

    def generate_dendrogram(
        self,
        mode: str,
        random_generator: np.random.Generator,
        *,
        alpha: float = 1.0,
        resolution: float = 1.0,
    ) -> None:
        current_graph_pos = self.ginfo_pos.graph_whole.copy()
        current_graph_neg = self.ginfo_neg.graph_whole.copy()
        partition_list: list[Partition] = []
        q = -1.0

        while True:
            self._move_nodes(random_generator, alpha=alpha, resolution=resolution)
            renumbered_partition = renumber(self.partition)
            partition_list.append(renumbered_partition)
            current_graph_pos = self._aggregate_nodes(renumbered_partition, current_graph_pos)
            current_graph_neg = self._aggregate_nodes(renumbered_partition, current_graph_neg)
            self.graph_whole = self.get_signed_graph(current_graph_pos, current_graph_neg, mode)

            self.partition = {node: i for i, node in enumerate(self.graph_whole.nodes())}
            self.ginfo_pos = GraphInfo(self.graph_whole, current_graph_pos, self.partition, self.weight_key)
            self.ginfo_neg = GraphInfo(self.graph_whole, current_graph_neg, self.partition, self.weight_key)
            new_q = self._modularity(self.partition, alpha=alpha, resolution=resolution)
            if new_q - q < MIN:
                break
            q = new_q

        self.dendrogram = partition_list[:]

    def best_partition(self, *, alpha: float = 1.0, resolution: float = 1.0, seed: int | None = None) -> Partition:
        if not isinstance(alpha, (int, float)):
            raise TypeError("Alpha value must be a number.")
        if not 0 <= float(alpha) <= 1:
            raise ValueError("Alpha value must be in the range [0, 1].")
        print(f"alpha: {alpha}, resolution: {resolution}")

        self.generate_dendrogram(self.mode, np.random.default_rng(seed=seed), alpha=float(alpha), resolution=resolution)
        partition = self.dendrogram[0].copy()
        for level in range(1, len(self.dendrogram)):
            for node, community in list(partition.items()):
                partition[node] = self.dendrogram[level][community]

        print(f"Final Modularity: {self._modularity(partition, alpha=float(alpha), resolution=resolution)}")
        print(f"Positive Modularity (Q+): {self.ginfo_pos.calc_modularity(partition, resolution=resolution)}")
        print(f"Negative Modularity (Q-): {self.ginfo_neg.calc_modularity(partition, resolution=resolution)}")
        return partition



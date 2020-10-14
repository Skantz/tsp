from typing import List, Tuple

import networkx as nx

from mst_graphs import Graph


class NetworkXGraph(Graph):
    """
    Implementation of `Graph` for local NetworkX graphs
    """

    graph: nx.Graph

    def __init__(self, graph: nx.Graph, size: int, max_weight: int, max_read: int):
        self.graph = graph

        super().__init__(size, max_weight, max_read)

    def _query_platform(self, node: int) -> List[Tuple[int, int]]:
        return [
            (neighbor, self.graph.edges[node, neighbor]['weight'])
            for neighbor in self.graph.neighbors(node)
        ]

# TODO: Probably easiest if we create Graph subclasses for each of the graphs we want to test

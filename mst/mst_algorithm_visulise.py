import random
import time
from queue import Queue
from typing import List

from mst_graphs import Graph, KattisGraph
import networkx as nx
import matplotlib.pyplot as plt

import matplotlib as mpl
from matplotlib import pyplot as plt
from pylab import rcParams
import networkx as nx

class MST:
    """<
    Class providing functions for the MST algorithm

    :ivar graph: `Graph` instance to query
    :ivar time_limit: Maximum time the algorithm is allowed to run
    :ivar start_time: Time the algorithm started (set at initialization)
    """

    graph: Graph
    time_limit: float
    start_time: float

    def __init__(self, graph: Graph, time_limit: float = 8.9):
        self.graph = graph
        self.time_limit = time_limit
        self.start_time = time.time()

    def approx_weight(self) -> float:
        """
        Approximate the weight of the minimum spanning forest
        """

        G = self.graph

        # Trivial cases
        if G.size <= 1 or G.max_weight == 0:
            return 0

        # Calculate component counts
        # These must be between 1 and the graph size, and can't be fractional
        components = [
            min(G.size, max(1, round(comp)))
            for comp in self.approx_components()
        ]

        weight = G.size + sum(components[:-1]) - components[-1] * G.max_weight

        # The weight must be between (size / 2) and ((size - 1) * max_weight))
        weight = max(G.size / 2, weight)
        weight = min(weight, (G.size - 1) * G.max_weight)

        return weight

    def approx_components(self) -> List[float]:
        """
        Approximate the number of components for each of the weights
        """

        G = self.graph

        # plt.figure(figsize=(18, 18))
        # nx.draw(G.graph, node_size=60)
        # plt.savefig("ba.png")
        # plt.show()

        color_lookup = {k: 0 for v, k in enumerate(sorted(set(G.graph.nodes())))}
        # {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'J': 6, 'K': 7, 'Z': 8}

        # nx.draw_networkx(G)
        # plt.show()
        low, *_, high = sorted(color_lookup.values())
        norm = mpl.colors.Normalize(vmin=low, vmax=high, clip=True)
        mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.coolwarm)

        rcParams['figure.figsize'] = 12, 7
        nx.draw(G.graph,
                nodelist=color_lookup,
                node_size=1000,
                node_color=[mapper.to_rgba(i)
                            for i in color_lookup.values()],
                with_labels=True)
        plt.show()


        i = 0

        # The number of components estimated for each weight restriction w <= 1, 2, ...
        betas = [0] * G.max_weight

        # Set to ensure the starting vertices are unique
        chosen_us = set()

        # Sample nodes until we hit the node limit
        while G.curr_read < G.max_read and time.time() - self.start_time < self.time_limit:
            i += 1

            while True:
                start = random.randint(0, G.size - 1)

                if start not in chosen_us:
                    chosen_us.add(start)
                    break

            # Randomly select the number of nodes sampled according to Pr[X] <= x = 1 - 1/x
            X = int(1.0 / (1.0 - random.random()))

            for weight in range(G.max_weight, 0, -1):
                nodes_queue = Queue()
                nodes_queue.put(start)
                queued_nodes = {start}

                while G.curr_read < G.max_read and time.time() - self.start_time < self.time_limit:
                    node = nodes_queue.get()
                    color_lookup[node] = 10

                    ####

                    low, *_, high = sorted(color_lookup.values())
                    norm = mpl.colors.Normalize(vmin=low, vmax=high, clip=True)
                    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.coolwarm)

                    rcParams['figure.figsize'] = 12, 7
                    nx.draw(G.graph,
                            nodelist=color_lookup,
                            node_size=1000,
                            node_color=[mapper.to_rgba(i)
                                        for i in color_lookup.values()],
                            with_labels=True)
                    plt.show()


                    ######ß


                    for neighbor in G.query(node):
                        if neighbor[0] not in queued_nodes and neighbor[1] <= weight:
                            nodes_queue.put(neighbor[0])
                            queued_nodes.add(neighbor[0])

                    if len(queued_nodes) > X:
                        break

                    if nodes_queue.empty():
                        betas[weight - 1] += 1
                        break

        if i <= 1:
            return [0 for _ in range(G.max_weight)]

        # If we find a component each time, i = betas[j] for all j.
        return [G.size / i * betas[j] for j in range(G.max_weight)]


if __name__ == "__main__":
    print("end", MST(KattisGraph()).approx_weight(), flush=True)

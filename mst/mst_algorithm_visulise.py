import random
import time
from queue import Queue
from typing import List

from mst.mst_graphs import Graph, KattisGraph
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
        print("approx_conponents is called\n")

        edge_list = [[], [], [], [], [], []]

        G = self.graph
        for i in range(5):
            G.graph.nodes[i]["x-coord"] = 10 * i
            G.graph.nodes[i]["y-coord"] = -1 * i * i + 4 * i + 4

        i = 0

        # The number of components estimated for each weight restriction w <= 1, 2, ...
        betas = [0] * G.max_weight

        # Set to ensure the starting vertices are unique
        chosen_us = set()

        # Sample nodes until we hit the node limit
        while G.curr_read < G.max_read and i < 25:
            i += 1
            print('Within time limit:')
            print(self.time_limit)
            print('i is:')
            print(i)

            while True:
                start = random.randint(0, G.size - 1)

                if start not in chosen_us:
                    chosen_us.add(start)
                    break

            # Randomly select the number of nodes sampled according to Pr[X] <= x = 1 - 1/x
            X = int(1.0 / (1.0 - random.random()))

            for weight in range(G.max_weight, 0, -1):

                ##### search for new weight, reset
                color_lookup = {k: 0 for v, k in enumerate(sorted(set(G.graph.nodes())))}
                labels_lookup = {}
                # {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'J': 6, 'K': 7, 'Z': 8}
                bfs_counter = 0
                #####

                print('############################# Searching weight ' + str(weight) + ' #############################')
                print('bfs counter is')
                print(bfs_counter)


                nodes_queue = Queue()
                nodes_queue.put(start)
                queued_nodes = {start}

                while G.curr_read < G.max_read and i < 25:
                    bfs_counter += 1

                    print('############################# Reading ' + str(G.curr_read) + ' #############################')
                    print('bfs counter is')
                    print(bfs_counter)

                    node = nodes_queue.get()

                    #### setting labels
                    labels_lookup[node] = str(bfs_counter)
                    #### end setting labels


                    # pos = nx.get_node_attributes(G.graph, 'pos')
                    # pos = nx.spring_layout(G.graph)  # positions for all nodes
                    # pos = dict()
                    # pos[0] = [0, 0]
                    # pos[1] = [2, 2]
                    # pos[2] = [2, 0]
                    # pos[3] = [0, 2]
                    # pos[4] = [1, 1]
                    atest = G.graph.nodes().data()
                    pos = {node[0]: (node[1]['x-coord'], node[1]['y-coord']) for node in G.graph.nodes(data=True)}

                    nx.draw(G.graph, pos=pos, node_size=20, alpha=0.3)
                    edge_labels = nx.get_edge_attributes(G.graph, 'weight')
                    nx.draw_networkx_edge_labels(G.graph, pos, edge_labels=edge_labels, alpha=0.3)
                    plt.show()
                    time.sleep(1)

                    node_to_query = node
                    for neighbor in G.query(node_to_query):
                        if neighbor[0] not in queued_nodes and neighbor[1] <= weight:

                            edge_list[weight].append([node_to_query, neighbor[0]])
                            node = neighbor[0]
                            print('Sampling node' + str(node))


                            bfs_counter += 1


                            #### draw
                            # node_colour = [mapper.to_rgba(i) for i in color_lookup.values()]
                            # labels = nx.get_edge_attributes(G.graph, 'weight')
                            # nx.draw_networkx_nodes(G.graph, pos, node_size=200, node_color=node_colour)
                            # nx.draw_networkx_edges(G.graph, pos, edgelist=G.graph.edges(data=True))
                            # nx.draw_networkx_labels(G.graph, pos, labels_lookup, font_size=20, font_color='white')
                            # nx.draw_networkx_edge_labels(G.graph, pos, edge_labels=labels)
                            # plt.show()
                            ######

                            nodes_queue.put(neighbor[0])
                            queued_nodes.add(neighbor[0])


                    if len(queued_nodes) > X:
                        break

                    if nodes_queue.empty():
                        betas[weight - 1] += 1
                        break
                print('############################# Finished BFS ############################# \n')
                print(weight)
                print(edge_list)

                # plot BFS for w == 5
                for i in range(len(edge_list[weight])):
                    edges_to_plot = edge_list[weight]
                    edges_to_plot = edges_to_plot[0:i+1]

                    g_odd_complete_min_edges = nx.Graph(edges_to_plot)
                    nx.draw(G.graph, pos=pos, node_size=20, alpha=0.05)
                    nx.draw(g_odd_complete_min_edges, pos=pos, node_size=20, edge_color='blue', node_color='red')

                    edge_labels = nx.get_edge_attributes(G.graph, 'weight')
                    nx.draw_networkx_edge_labels(G.graph, pos, edge_labels=edge_labels, alpha=0.3)
                    print('w ==' + str(weight))
                    plt.text(0, 7.8, "weight == "+str(weight))
                    plt.show()
                    time.sleep(1)
                # Finish plot BFS for w == 5

        print('time up')
        print(time.time() - self.start_time)
        if i <= 1:
            return [0 for _ in range(G.max_weight)]

        # If we find a component each time, i = betas[j] for all j.
        return [G.size / i * betas[j] for j in range(G.max_weight)]


if __name__ == "__main__":
    print("end", MST(KattisGraph()).approx_weight(), flush=True)

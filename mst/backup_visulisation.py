# First networkx library is imported
# along with matplotlib
import networkx as nx
import matplotlib.pyplot as plt

import matplotlib as mpl
from matplotlib import pyplot as plt
from pylab import rcParams
import networkx as nx

# Defining a Class
class GraphVisualization:

    def __init__(self):
        # visual is a list which stores all
        # the set of edges that constitutes a
        # graph
        self.visual = []

        # addEdge function inputs the vertices of an

    # edge and appends it to the visual list
    def addEdge(self, a, b):
        temp = [a, b]
        self.visual.append(temp)

        # In visualize function G is an object of

    # class Graph given by networkx G.add_edges_from(visual)
    # creates a graph with a given list
    # nx.draw_networkx(G) - plots the graph
    # plt.show() - displays the graph
    def visualize(self):
        G = nx.Graph()
        G.add_edges_from(self.visual)
        # color_lookup = {0:0, 1:0, 2:0, 3:0, 4:0, 5:1}

        color_lookup = {k: 0 for v, k in enumerate(sorted(set(G.nodes())))}
        # {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'J': 6, 'K': 7, 'Z': 8}

        # nx.draw_networkx(G)
        # plt.show()
        low, *_, high = sorted(color_lookup.values())
        norm = mpl.colors.Normalize(vmin=low, vmax=high, clip=True)
        mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.coolwarm)

        rcParams['figure.figsize'] = 12, 7
        nx.draw(G,
                nodelist=color_lookup,
                node_size=1000,
                node_color=[mapper.to_rgba(i)
                            for i in color_lookup.values()],
                with_labels=True)
        plt.show()


    # Driver code


G = GraphVisualization()
G.addEdge(0, 2)
G.addEdge(1, 2)
G.addEdge(1, 3)
G.addEdge(5, 3)
G.addEdge(3, 4)
G.addEdge(1, 0)
G.visualize()
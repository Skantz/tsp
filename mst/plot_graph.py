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
import pandas as pd
import copy

def normalize(edge):
    (n1, n2, dist) = edge
    if n1 > n2: # use a custom compare function if desired
        n1, n2 = n2, n1
    return ((n1, n2, dist))


edge_list = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
g = nx.nx.Graph()
g.add_weighted_edges_from(
    [(0, 1, 3), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 2, 1), (1, 3, 2), (1, 4, 3), (2, 3, 2), (2, 4, 2), (3, 4, 2)])

g_old = copy.deepcopy(g)
# pos = nx.spring_layout(g)
# print(pos)
pos = {0: ([0.2747637, -1.]), 1: ([-0.52459149, -0.48841505]), 2: ([0.24988574, 0.96947934]),
       3: ([0.66764649, 0.08829556]), 4: ([-0.66770444, 0.43064015])}

labels_lookup = {}

### draw
labels = nx.get_edge_attributes(g, 'weight')
nx.draw_networkx_nodes(g, pos, node_size=200, node_color='gray')
nx.draw_networkx_edges(g, pos, edgelist=g.edges(data=True))
# nx.draw_networkx_labels(g, pos, labels_lookup, font_size=20, font_color='white')
nx.draw_networkx_edge_labels(g, pos, edge_labels=labels)
plt.show()
######

# bug in plotting
# nx.draw_networkx_edges(g, pos, edgelist=[e for e in g.edges(data=True) if e[2]['weight'] < 3])
# nx.draw_networkx_labels(g, pos, labels_lookup, font_size=20, font_color='white')
# nx.draw_networkx_edge_labels(g, pos, edge_labels={key: labels[key] for key in labels if g.edges[key]['weight'] < 3})
# plt.show()


# w = 3
### draw whole graph
labels = nx.get_edge_attributes(g, 'weight')
nx.draw_networkx_nodes(g, pos, node_size=200, node_color='gray')
nx.draw_networkx_edges(g, pos, edgelist=g.edges(data=True), alpha=0.2)
# nx.draw_networkx_labels(g, pos, labels_lookup, font_size=20, font_color='white')
nx.draw_networkx_edge_labels(g, pos, edge_labels=labels, alpha=0.2)
######


remove_list = []
for n, nbrs in g.adjacency():
    for nbr, eattr in nbrs.items():
        data = eattr['weight']
        if data >= 3:
            print('(%d, %d, %0.3f)' % (n, nbr, data))
            remove_list.append((n, nbr, data))
print('***********************************')

unique_edges = set(map(normalize, remove_list))

for e in unique_edges:
    g.remove_edge(e[0], e[1])


labels = nx.get_edge_attributes(g, 'weight')
nx.draw_networkx_nodes(g, pos, node_size=200, node_color='gray')
nx.draw_networkx_edges(g, pos, edgelist=g.edges(data=True))
# nx.draw_networkx_labels(g, pos, labels_lookup, font_size=20, font_color='white')
nx.draw_networkx_edge_labels(g, pos, edge_labels=labels)
plt.show()


bfs_path_edges = [[0, 2], [0, 3], [0, 4], [4, 1], [4, 2], [4, 3]]
g_odd_complete_min_edges = nx.Graph(bfs_path_edges)
labels = nx.get_edge_attributes(g_odd_complete_min_edges, 'weight')
nx.draw_networkx_nodes(g_odd_complete_min_edges, pos, node_size=200, node_color='red')

# draw old graph
labels = nx.get_edge_attributes(g, 'weight')
nx.draw_networkx_nodes(g, pos, node_size=200, node_color='gray', alpha=0.2)
nx.draw_networkx_edges(g, pos, edgelist=g.edges(data=True), alpha=0.2)
# nx.draw_networkx_labels(g, pos, labels_lookup, font_size=20, font_color='white')
nx.draw_networkx_edge_labels(g, pos, edge_labels=labels)
# finish old graph

nx.draw_networkx_edges(g_odd_complete_min_edges, pos, edgelist=g_odd_complete_min_edges.edges(data=True), edge_color='blue')
# nx.draw_networkx_labels(g, pos, labels_lookup, font_size=20, font_color='white')
nx.draw_networkx_edge_labels(g, pos, edge_labels=labels)
plt.show()


# w = 2
### draw whole graph
labels = nx.get_edge_attributes(g, 'weight')
nx.draw_networkx_nodes(g, pos, node_size=200, node_color='gray')
nx.draw_networkx_edges(g, pos, edgelist=g.edges(data=True), alpha=0.2)
# nx.draw_networkx_labels(g, pos, labels_lookup, font_size=20, font_color='white')
nx.draw_networkx_edge_labels(g, pos, edge_labels=labels, alpha=0.2)
######

# draw subgraph
remove_list = []
for n, nbrs in g.adjacency():
    for nbr, eattr in nbrs.items():
        data = eattr['weight']
        if data >= 2:
            print('(%d, %d, %0.3f)' % (n, nbr, data))
            remove_list.append((n, nbr, data))
print('***********************************')

unique_edges = set(map(normalize, remove_list))

for e in unique_edges:
    g.remove_edge(e[0], e[1])


labels = nx.get_edge_attributes(g, 'weight')
nx.draw_networkx_nodes(g, pos, node_size=200, node_color='gray')
nx.draw_networkx_edges(g, pos, edgelist=g.edges(data=True))
# nx.draw_networkx_labels(g, pos, labels_lookup, font_size=20, font_color='white')
nx.draw_networkx_edge_labels(g, pos, edge_labels=labels)
plt.show()

time.sleep(1)

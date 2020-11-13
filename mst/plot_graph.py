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

edge_list = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
g = nx.nx.Graph()
g.add_weighted_edges_from(
    [(0, 1, 3), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 2, 1), (1, 3, 2), (1, 4, 3), (2, 3, 2), (2, 4, 2), (3, 4, 2)])

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


time.sleep(1)

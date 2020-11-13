import matplotlib.pyplot as plt
import numpy as np
from mst.mst_algorithm_visulise import MST
from mst.mst_graphs_local import NetworkXGraph
import networkx as nx
import random
import time
from enum import Enum


# Constants

class GraphType(Enum):
    RANDOM = 'ramdom_graph'
    EXPANDER = 'expander_graph'
    COMPLETE = 'complete_graph'
    CYCLE = 'cycle_graph'
    VISUAL = 'visual_graph'


class WeightDistribution(Enum):
    RANDOM = 'random_weights'
    EXPONENTIAL = 'exponantial_weights'
    NORMAL = 'normal_weights'
    VISUAL = 'visual_weights'


default_size = 100
default_degree = 2
default_max_weight = 100
default_max_read = 30
default_graph_type = GraphType.RANDOM
default_weight_distribution = WeightDistribution.RANDOM


def print_constants(custom_size, degree, max_read, max_weight, graph_type, weight_distribution):
    output = "Size: " + str(custom_size) + "\n" + "Degree: " + str(degree) + "\n" + "Max Weight: " + str(
        max_weight) + "\n" + "Max Read: " + str(max_read) + "\n" + "Graph Type: " + str(
        graph_type) + "\n" + "Weight Distribution: " + str(weight_distribution) + "\n"
    print(output)


# %% Graph generation functions

# Random NetworkX graph

def random_graph(degree, size):
    graph = nx.generators.random_regular_graph(degree, size)
    return graph


# Expanders NetworkX Graph

def expander_graph(size):
    graph = nx.generators.paley_graph(size)
    return graph


# Complete graph

def complete_graph(size):
    graph = nx.generators.complete_graph(size)
    return graph


# Cycle graph

def cycle_graph(size):
    graph = nx.generators.cycle_graph(size)
    return graph

# Visual graph

def visual_graph(size):
    graph = nx.generators.complete_graph(size)
    return graph

# Assign edge weights

def assign_random_edge_weights(graph, max_weight):
    for (u, v) in graph.edges():
        graph.edges[u, v]['weight'] = random.randint(1, max_weight)
    return graph

def assign_exponential_edge_weights(graph,custom_scale=1):
    for (u, v) in graph.edges():
        graph.edges[u, v]['weight'] = round(np.random.exponential(scale=custom_scale)+1)
    return graph

def assign_normal_edge_weights(graph, max_weight):
    for (u, v) in graph.edges():
        graph.edges[u, v]['weight'] = round(np.random.normal(loc=(max_weight/2)))
    return graph

def assign_visual_edge_weights(graph, max_weight):
    for (u, v) in graph.edges():
        graph.edges[u, v]['weight'] = random.randint(0, 9) % 3 + 1
    return graph

# Running against nx.kruskall

def kruskall_comparison(custom_size=default_size, degree=default_degree, max_read=default_max_read,
                        max_weight=default_max_weight, graph_type=default_graph_type,
                        weight_distribution=default_weight_distribution):
    if max_read >= custom_size:
        raise Exception("Please reduce default_max_read or increase node size!")

    print("Create graph with the following constants:")
    print_constants(custom_size, degree, max_read, max_weight, graph_type, weight_distribution)

    # Chose here different graph options from above
    graph = None

    if graph_type is GraphType.RANDOM:
        graph = random_graph(degree, custom_size)
    elif graph_type is GraphType.EXPANDER:
        graph = expander_graph(custom_size)
    elif graph_type is GraphType.COMPLETE:
        graph = complete_graph(custom_size)
    elif graph_type is GraphType.CYCLE:
        graph = cycle_graph(custom_size)
    elif graph_type is GraphType.VISUAL:
        graph = visual_graph(custom_size)

    # Chose here different weight sampling options from above
    if weight_distribution is WeightDistribution.RANDOM:
        graph = assign_random_edge_weights(graph, max_weight)
    elif weight_distribution is WeightDistribution.EXPONENTIAL:
        graph = assign_exponential_edge_weights(graph)
    elif weight_distribution is WeightDistribution.NORMAL:
        graph = assign_normal_edge_weights(graph, max_weight)
    elif weight_distribution is WeightDistribution.VISUAL:
        graph = assign_visual_edge_weights(graph, max_weight)

    # Our algorithm
    our_start_time = time.time()
    approx = MST(
        NetworkXGraph(graph, custom_size, max_weight, max_read)
    ).approx_weight()
    our_runtime = time.time() - our_start_time
    # print("Our Runtime: "+str(our_runtime))
    # print("Our MST weight approx. for random graph:", approx,"\n")

    # NetworkX Kruskall Algorithm
    kruskall_start_time = time.time()
    mst = nx.minimum_spanning_tree(graph)
    kruskall_time = time.time() - kruskall_start_time
    # print("NetworkX Kruskall runtime: "+str(kruskall_time))

    gt_w = sum(graph.edges[e]["weight"] for e in mst.edges)
    # print("ground truth weight", gt_w)
    print(gt_w)
    print(approx)

    ratio = max(gt_w, approx) / min(gt_w, approx)
    err = 100 * (ratio - 1)
    # print("error %", err)

    results = {
        "our_time": our_runtime,
        "kuskrall_time": kruskall_time,
        "groundtruth": gt_w,
        "approximation": approx,
        "error": err
    }

    print("Results", results, "\n")
    return results

result_single_custom = kruskall_comparison(custom_size=5, degree=2, max_read=3, max_weight=3, graph_type=GraphType.VISUAL,weight_distribution=WeightDistribution.VISUAL)









# end
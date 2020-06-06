#
# Example usage of tsp_local on a small subset of U.S. cities, more can be
# found in the TSPLib, fri26.tsp.
#
from tsp_local.test import matrix
from tsp_local.base import TSP
from tsp_local.kopt import KOpt

names = [
    "New York",
    "Los Angeles",
    "Chicago",
    "Minneapolis",
    "Denver",
    "Dallas",
    "Seattle",
    "Boston",
    "San Francisco",
    "St. Louis",
    "Houston",
    "Phoenix",
    "Salt Lake City"
]
# Load the distances
TSP.setEdges(matrix)
# Make an instance with all nodes
lk = KOpt(range(len(matrix)))
path, cost = lk.optimise()
print("Best path has cost: {}".format(cost))
print([names[i] for i in path])

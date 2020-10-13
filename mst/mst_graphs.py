from abc import ABC, abstractmethod
from typing import Dict, List, Tuple


class Graph(ABC):
    """
    Abstract class that provides access to graphs to the MST algorithm

    :ivar size: Number of nodes in the graph
    :ivar max_weight: Maximum edge weight in the graph
    :ivar curr_read: Current number of edges read from the graph
    :ivar max_read: Maximum number of edges to be read from the graph
    :ivar memos: Memorized values for queried nodes
    """

    size: int
    max_weight: int
    curr_read: int
    max_read: int
    memos: Dict[int, List[Tuple[int, int]]]

    def __init__(self, size: int, max_weight: int, max_read: int):
        self.size = size
        self.max_weight = max_weight
        self.curr_read = 0
        self.max_read = max_read

        self.memos = {}

    def query(self, node: int) -> List[Tuple[int, int]]:
        """
        Query a node for its incident edges

        :param node: Index of the requested node
        :return: List of incident edges with for each the neighbor node and weight
        """

        if node not in self.memos:
            self.memos[node] = self._query_platform(node)
            self.curr_read += 1

        return self.memos[node]

    @abstractmethod
    def _query_platform(self, node: int) -> List[Tuple[int, int]]:
        pass


class KattisGraph(Graph):
    """
    Implementation of `Graph` for Kattis-provided graphs
    """

    def __init__(self):
        size = int(input())
        max_weight = int(input())
        max_read = int(input())

        super().__init__(size, max_weight, max_read)

    def _query_platform(self, node: int) -> List[Tuple[int, int]]:
        print(node, flush=True)

        result = input().split()

        return [
            (int(result[i]), int(result[i + 1]))
            for i in range(1, len(result), 2)
        ]

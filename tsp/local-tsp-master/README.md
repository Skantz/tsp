# TSP Local

Python implementation of different local heuristics for the TSP. Was developped
as part of a larger project so it might look quirky.

The basic idea is to have a parent class holding the data for a TSP instance
and then create _scenarios_ based on this instance and solve them. Obviously
works as a one-off.

## Usage

A typical workflow is:

 0. Assign the distance matrix with `setEdges()`.
 1. Create the scenario with `__init__()`.
 2. At every iteration, use `update()` to set the `heuristic_*` fields.
 3. Then improve the path using `optimise()`, calling the local improvement.

```python
# Load the distances
TSP.setEdges(matrix)
# Make an instance with all nodes
lk = KOpt(range(len(matrix)))
path, cost = lk.optimise()
```

## Installation

This is a standalone library, better included directly in projects. Works with
Python 2 and 3 hopefully.

### Optional

 * [Doctest][doctest] for running tests.
 * [Pytest][pytest] for test integration.
 * [Yapf][yapf] for code formatting.

### Contributing

This project has a few additional requirements which should be enforced inside
`setup.cfg`.

## TODO

 * Improve LKH, some of Helsgaun work is in there but more can be added.
 * Verify that both Python 2 and 3 work -- currently 3.5 okay for sure.

[doctest]: https://docs.python.org/2/library/doctest.html
[pytest]: https://docs.pytest.org/en/latest/
[yapf]: https://github.com/google/yapf
[tsplib]: https://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/
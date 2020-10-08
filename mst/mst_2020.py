
import random
import sys
import time
from queue import Queue
from math import log10, log2

#If we want to read local graphs (refer to the ipynb nobtebook)
DEBUG = 0

#Not signifcant
RANDOM_SAMPLE_FACTOR = 0.5
SOME_BIG_CONSTANT    = 1500

#Memoize to save reads
memos = {}
#Unused
sg_memos = {}

#Tracking read node counts from all functions
global read_nodes
read_nodes = 0

# See notebook for local testing
def query_node_local(n):
    """ Read a node locally from a NetworkX graph G """
    nbors = G.neighbors(n)
    return [(n2, G.edges[n, n2]["weight"]) for n2 in nbors]


#Return memoization if exists, otherwise ask for a new node
def query_node(n):
    """Return memoization of neighbors if exists, otherwise query node from stdio"""

    if DEBUG:
        return query_node_local(n)

    if n in memos:
        return memos[n]

    else:
        print(n)
        sys.stdout.flush()
        inp = input().strip("\n").split()[1:]
        memos[n] = []
        if len(inp) < 2:
            assert(len(inp) == 0)
            return memos[n]

        for i in range(0, len(inp), 2):

            #A neigbor list of edges : tuples (node, weight)
            memos[n].append((int(inp[i]), int(inp[i + 1])))

        global read_nodes        
        read_nodes += 1
        return memos[n]



def approx_cc_simple(n, gi, eps, max_w, d_bar, max_nodes, time_limit=None, start_time=None, read_new_verts=True, sample_size=0.5):
    """ Approximate a minimum spanning forest """

    #If eps > 1, still sample some constant or fraction

    betas = 0
    
    i = 0

    #The numbber of components estimated for each weight restriction w <= 1, 2, ...
    betas = [0 for _ in range(max_w)]

    #read_nodes = 0

    sampled = False

    # Sample nodes until we hit the node limit
    while True:
        if (time.time() - start_time) > time_limit:
            break

        i += 1

        # Unused. We always read new vertices.
        u = random.randint(0, n - 1)

        if read_nodes >= max_nodes:
            i -= 1
            break

        #Randomly select the number of nodes sampled according to Pr[X] <= x = 1 - 1/x
        X = int(1/random.random())

        for w in range(max_w, 0, -1):
            visited_nodes = set()
            visited_nodes.add(u)

            u_nbors_queue = Queue()
            u_nbors = [v for v in query_node(u) if v[0] != u and v[1] <= w]
            u_nbors = list(set(u_nbors))

            if read_nodes >= max_nodes:
                #i -= 1    #Off-by-ones
                break

            du = len(u_nbors)
            seen_edges = set([(u, v[0]) for v in u_nbors])
            for nb in u_nbors:
                u_nbors_queue.put(nb)
            
            if DEBUG:
                for v in u_nbors:
                    assert(u_nbors.count(v) == 1)

            if u_nbors == []:
                #betas[i] = 1
                betas[w - 1] += 1
                continue

            j = -1
            heads = 0
            
            #Constant lookup for processed nodes
            added_to_queue = set([u])
           
            #bfs
            while not u_nbors_queue.empty():
                if (time.time() - start_time) > time_limit:
                    break

                #v = u_nbors.pop(0)
                v = u_nbors_queue.get()
                j += 1
                visited_nodes.add(v[0])
                v_nbors = query_node(v[0])
                #v_nbors = sorted(v_nbors, key=lambda x: x[0] in memos, reverse=True)
                if read_nodes >= max_nodes:
                    #i -= 1
                    break
                added = False
                for t in v_nbors:
                    seen_edges.add((v[0], t[0]))
                    if t[0] not in visited_nodes and t[0] not in added_to_queue and t[1] <= w:
                        u_nbors_queue.put(t)
                        added_to_queue.add(t[0])


                #One of these conditions implies the others
                if j > X or len(visited_nodes) > X or len(added_to_queue) > X:
                    break            

                if u_nbors_queue.empty():
                    betas[w - 1] += 2**heads 

                    #1/2 chance to keep going. This might improve accuracy on dense graphs
                    toss = random.randint(0, 1)
                    if toss == 0:
                        break
                    else:
                        X *= 2
                        heads += 1

    if i <= 1:
        return [0 for _ in range(max_w)]

    #If we find a component each time, i = betas[j] for all j.
    return [n / i * betas[j] for j in range(max_w)]


def approx_mst(n, eps, max_w, max_nodes):

    d_bar = 0
    c_bars = []

    max_nodes = max_nodes  * 1

    #Several of these parameters are not used.
    c_bars = approx_cc_simple(n, max_w, eps, max_w, d_bar, max_nodes, time_limit=8.8, start_time=time.time(), read_new_verts=True)

    #We can't have fractional components
    c_bars = [min(n, max(1, round(bar))) for bar in c_bars]
    #We have at least one component without weight restriction
    c_bars[-1] = max(c_bars[-1], 1)

    if DEBUG == 1:
        print("est components", c_bars[-1])
        print("actual components", c)
        print("raw score", n  + sum(c_bars[1:]))

    return n + sum(c_bars[:-1]) - c_bars[-1] * max_w 


if not DEBUG:
    n = int(input()) # -> 2^32
    eps = 0.1
    max_w = int(input()) # 1 -> 50
    max_nodes = int(input()) - 1


#Trivial cases
if n < 2:
    print("end " + str(0))
if max_w < 1:
    print("end " + str(0))


mst = approx_mst(n, eps, max_w, max_nodes)

#A minimum spanning tree can't have less or more weight than this
mst = max(n/2, mst)
mst = min((n - 1) * max_w, mst)

print("end " + str(mst))
sys.stdout.flush()

#Debug prints for local tests
if DEBUG == 1:
    print("gt end", mst_gt_w)

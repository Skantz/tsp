
import random
import sys
import time
from queue import Queue
from math import log10, log2

DEBUG = 0

RANDOM_SAMPLE_FACTOR = 0.5
SOME_BIG_CONSTANT    = 1500

#Vertex to [(neighbor, edge weight), ... ]
memos = {}
sg_memos = {}


def query_node_local(n):
    nbors = G.neighbors(n)
    return [(n2, G.edges[n, n2]["weight"]) for n2 in nbors]


def query_node_in_subgraph_n(n, w):
    if (n, w) not in memos:
        u_nbors = [v for v in query_node(n) if v[0] != n and v[1] <= w]
        u_nbors = list(set(u_nbors))
        memos[(n, w)] = u_nbors
        return memos[(n, w)]
    else:
        return memos[(n, w)]


def query_node(n):

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

            #print("i", i, "len inp", len(inp), "inp", inp)
            memos[n].append((int(inp[i]), int(inp[i + 1])))
        
        return memos[n]



def approx_cc_simple(n, gi, eps, max_w, d_bar, max_nodes, time_limit=None, start_time=None, read_new_verts=True, sample_size=0.5):

    #If eps > 1, still sample some constant or fraction
    r = 5*min(1000, log10(n)) + int(RANDOM_SAMPLE_FACTOR/eps**2) #RANDOM_SAMPLE_FACTOR)

    betas = 0
    
    i = 0
    betas = [0 for _ in range(max_w)]

    read_nodes = 0

    sampled = False

    while i < r:
        if (time.time() - start_time) > time_limit:
            break

        #if i > 0 and (i % r*sample_size) == 0 and sampled == False:
        if i == r - 1:
            sampled = True
            checkpoint_score = sum([n / i * betas[j] for j in range(max_w)])
            if checkpoint_score < n/2 or checkpoint_score > max_w * n:
                r *= 2
                continue
            #r *= max(1, log2(betas[-1]))


        i += 1

        if read_new_verts:
            u = random.randint(0, n - 1)
        else:
            u = random.choice(list(memos.keys()))

        read_nodes += 1
        if read_nodes >= max_nodes:
            break

        X = int(1/random.random())

        #betas[i] = 0
        for w in range(max_w, 0, -1):
            visited_nodes = set()
            visited_nodes.add(u)

            u_nbors_queue = Queue()
            u_nbors = [v for v in query_node(u) if v[0] != u and v[1] <= w]
            u_nbors = list(set(u_nbors))

            read_nodes += 1
            if read_nodes >= max_nodes:
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
            
            while not u_nbors_queue.empty():
                if (time.time() - start_time) > time_limit:
                    break

                #v = u_nbors.pop(0)
                v = u_nbors_queue.get()
                j += 1
                visited_nodes.add(v[0])
                v_nbors = query_node(v[0])
                #v_nbors = sorted(v_nbors, key=lambda x: x[0] in memos, reverse=True)
                added = False
                for t in v_nbors:
                    seen_edges.add((v[0], t[0]))
                    if t[0] not in visited_nodes and t[0] not in added_to_queue and t[1] <= w:
                        u_nbors_queue.put(t)
                        added_to_queue.add(t[0])



                if j > X or len(visited_nodes) > X or len(added_to_queue) > X:
                    break            

                if u_nbors_queue.empty():
                    betas[w - 1] += 2**heads 

                    toss = random.randint(0, 1)
                    if toss == 0:
                        break
                    else:
                        X *= 2
                        heads += 1

    if i <= 1:
        return [0 for _ in range(max_w)]

    return [n / i * betas[j] for j in range(max_w)]


def approx_mst(n, eps, max_w, max_nodes):

    d_bar = 0
    c_bars = []

    #eps = eps / (1 + 0.5 * max_w)
    eps = eps / (1 + 0.5 * max_w)
    max_nodes = max_nodes  * 10

    if max_w > 1:
        #c_bars = approx_cc_simple(n, max_w, eps/2.5, max_w, d_bar, time_limit=2.9, start_time=time.time(), read_new_verts=True)
        c_bars = approx_cc_simple(n, max_w, eps, max_w, d_bar, max_nodes, time_limit=8.8, start_time=time.time(), read_new_verts=True)
    else:
        c_bars = approx_cc_simple(n, max_w, eps/10, max_w, d_bar, max_nodes, time_limit=8.8, start_time=time.time(), read_new_verts=True)

    c_bars = [min(n, max(1, round(bar))) for bar in c_bars]
    c_bars[-1] = max(c_bars[-1], 1)

    if DEBUG == 1:
        print("est components", c_bars[-1])
        print("actual components", c)
        print("raw score", n  + sum(c_bars[1:]))

    return n + sum(c_bars[:-1]) - c_bars[-1] * max_w 


if not DEBUG:
    n = int(input())
    #eps = float(input()) - 1
    eps = 0.1
    max_w = int(input())
    max_nodes = int(input()) - 1

#All groups have some eps > 1 (first 6 problems)
#Probably eps is between 0 and 0.5


if n < 2:
    print("end " + str(0))

if max_w < 1:
    print("end " + str(0))


#Group 6 has max_w > 50
#if max_w > 50:
#    1/0


mst = approx_mst(n, eps, max_w, max_nodes)
mst = max(n/2, mst)
mst = min((n - 1) * max_w, mst)

print("end " + str(mst))
sys.stdout.flush()


if DEBUG == 1:
    print("gt end", mst_gt_w)

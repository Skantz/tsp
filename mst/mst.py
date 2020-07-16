import random
import sys

DEBUG = 0

RANDOM_SAMPLE_FACTOR = 200
SOME_BIG_CONSTANT    = 1000

#Vertex to [(neighbor, edge weight), ... ]
memos = {}

def query_node_local(n):
    nbors = G.neighbors(n)
    return [(n2, G.edges[n, n2]["weight"]) for n2 in nbors]

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
            return memos[n]

        for i in range(0, len(inp), 2):
            #print("i", i, "len inp", len(inp), "inp", inp)
            memos[n].append((int(inp[i]), int(inp[i + 1])))
        
        return memos[n]

def approx_avg_degree(eps, C, n):
    import random
    n_ = 50 + int(C/eps)
    n = min(n_, 1+ n//1000)
    max_deg = 0

    for i in range(n_):
        node = random.randint(0, n - 1)
        out = query_node(node)

        deg = len(out)
        max_deg = max(max_deg, deg)


    return max_deg


def approx_cc(n, gi, eps, max_w, d_bar):
    import random
    import sys

    #if DEBUG:
    #    print("gi", gi, "eps", eps, "max_w", max_w, "d bar", d_bar)

    r = int(1/eps) * RANDOM_SAMPLE_FACTOR
    #Does not check dupes
    vs = [random.randint(0, n - 1) for _ in range(r)]
    for v in vs:
        for _ in range(vs.count(v) - 1):
            vs.remove(v)
    
    r = len(vs)
    betas = [0 for _ in range(int(r))]

    for i, v in enumerate(vs):
        nb_w_tups = query_node(v)
        nb_w_tups = [t for t in nb_w_tups if t[1] <= gi]
        new_n = [t for t in nb_w_tups]
        n_visit = 0
        max_degree_seen = 0
        nodes_visited = set()
        edges_visited = set()

        if len(nb_w_tups) == 0:
            singleton = True
        else:
            singleton = False
        #singleton = True if len(nb_w_tups) == 0 else False
        heads = 0
        dui = len(nb_w_tups)
        walk_limit = max_w
        
        e_visit = 0

        for j, nb_w in enumerate(nb_w_tups):
            "BFS"
            #bfs
            nodes_visited.add(nb_w[0])
            for e in [(v, tup[0]) for tup in nb_w_tups]:
                edges_visited.add(e)

            toss = random.randint(0, 1)
            heads += 1
            max_degree_seen = max(max_degree_seen, len(memos[nb_w[0]]))
            #print("toss", toss, "n_visit", n_visit, "walk limit", walk_limit, "max deg seen", max_degree_seen, "d bar", d_bar)
            if toss and len(edges_visited) <= walk_limit and max_degree_seen <= d_bar:
                for tup in nb_w_tups[1]:
                    new_n = query_node(tup[0])
                    nb_w_tups += [nn for nn in new_n if nn[1] <= gi and nn[0] not in nodes_visited and nn not in nb_w_tups]
                    #walk_limit = walk_limit*2
                #max_degree_seen = max(max_degree_seen, dui)
                n_visit+=1
                walk_limit *= 2
            else:
                break

        if singleton:
            betas[i] = 2
            continue
        elif len(nb_w_tups) == 0:
            #Finished bfs
            betas[i] = dui * 2**(heads)
   
    return n / (2 * r) * sum(betas)


def approx_cc_simple(n, gi, eps, max_w, d_bar):
    import random
    import sys

    #r = min(max(1, int(1/eps) * min(RANDOM_SAMPLE_FACTOR, 1 + int(n/RANDOM_SAMPLE_FACTOR))), n)
    r = 50 + int(1/eps * RANDOM_SAMPLE_FACTOR)
    r = min(1 + n//1000, r)
    #Does not check dupes
    vs = [random.randint(0, n - 1) for _ in range(r)]

    betas = [1 for _ in range(int(r))]


    #BFS
    queue = []
    
    for i, u in enumerate(vs):
        X = int( 1/random.random())

        betas[i] = 0
        visited_nodes = set()
        visited_nodes.add(u)

        u_nbors = [v for v in query_node(u) if v[0] != u and v[1] <= gi]

        for v in u_nbors:
            for _ in range(u_nbors.count(v) - 1):
                u_nbors.remove(v)

        if DEBUG:
            for v in u_nbors:
                assert(u_nbors.count(v) == 1)

        if u_nbors == []:
            betas[i] = 1
            continue

        j = 0
        for v in u_nbors:
            j += 1
            visited_nodes.add(v[0])
            v_nbors = query_node(v[0])
            added = False
            for t in v_nbors:
                if t[0] not in visited_nodes and t[0] not in u_nbors and t[1] <= gi:
                    u_nbors.append(t)
                    added = True

            if j > X  or len(visited_nodes) > X :
                betas[i] = 0
                break
            
            if u_nbors[j:] == [] or not added:
                betas[i] = 1
                break


    return n / r * sum(betas)


def approx_mst(n, eps, max_w):

    d_bar = approx_avg_degree(eps, SOME_BIG_CONSTANT, n)

    c_bars = []
    for i in range(1, max_w):

        cc_bar = approx_cc_simple(n, i , eps, 4*max_w/eps, d_bar)   #G^i : edges with weight at most i
        c_bars.append(c
    if DEBUG == 1:
        print("predicted n components", all_)
        print("actual    n components", c)c_bar)

    all_ = approx_cc_simple(n, max_w, eps/10, 4*max_w/eps, d_bar)
    if DEBUG == 1:
        print("estimated n components", all_)
    return n - max_w + sum(c_bars) - all_*max_w #, approx_cc_simple(n, max_w, eps, 4*max_w/eps, d_bar)

#Doesn't work outside ipython
if not DEBUG:
    n = int(input())
    eps = float(input())
    max_w = int(input())

mst = approx_mst(n, eps, max_w)

#assert(mst >= n/2)

print("end " + str(mst))
sys.stdout.flush()

if DEBUG == 1:
    print("gt end", mst_gt_w)

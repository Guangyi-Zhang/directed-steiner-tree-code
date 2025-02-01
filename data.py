import numpy as np
import networkx as nx

np.random.seed(42)


weighted = True
nT = 10
nNeighbor = 10

#for n in [1e3,1e4,1e5,1e6]:
#for n in [1e7,]:
for n in [1e8,]:
    n = int(n)
    p = nNeighbor/n
    g = nx.fast_gnp_random_graph(n, p, directed=True, seed=42)
    print(f"{n} done")

    edges = []
    edgeweights = []
    for v, nbrs in g.adj.items(): # n starts with 0
        edges.extend([(v, nbr) for nbr, eattr in nbrs.items()]) # eattr['weight']
        if weighted:
            edgeweights.extend(np.random.rand(len(nbrs)))
        else:
            edgeweights.extend([1] * len(nbrs))

    sp = str(p).replace('.','')
    with open(f"datasets/random_graph_w{weighted}_nghb{nNeighbor}_n{n}.csv", 'w') as fout:
        for e,w in zip(edges, edgeweights):
            fout.write(f"{e[0]},{e[1]},{w}")
            fout.write('\n')

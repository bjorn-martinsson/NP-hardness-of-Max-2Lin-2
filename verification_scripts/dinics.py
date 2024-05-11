# Implementation of Dinics max flow algorithm
# by pajenegod (@bjorn-martinsson on github)

class flow_graph:
    def __init__(self, add_edge, Dinic, calc_flow):
        self.add_edge = add_edge
        self.solve = Dinic
        self.calc_flow = calc_flow

def init():
    V = []
    C = []
    F = []
    def add_edge(u, v, cap, rcap = 0):
        V.append(v)
        V.append(u)
        C.append(cap)
        C.append(rcap)

    def Dinic(s, t):
        n = max(V) + 1
        m = len(V)
        coupl = [[] for _ in range(n)]
        for i in range(m):
            coupl[V[i]].append(i ^ 1)
        F[:] = C

        P = [-1]*n
        total_flow = 0
        while True:
            order = [0]*n
            bfs = [s]
            order[s] = 1
            for node in bfs:
                for eind in coupl[node]:
                    if not order[V[eind]] and F[eind]:
                        bfs.append(V[eind])
                        order[V[eind]] = len(bfs)
                        if V[eind] == t:
                            break
                else:
                    continue
                break
            else:
                break
            
            edge_counter = [len(c) for c in coupl]
            node = s
            while True:
                while edge_counter[node]:
                    eind = coupl[node][edge_counter[node] - 1]
                    if order[V[eind]] > order[node] and F[eind]:
                        node = V[eind]
                        P[node] = eind ^ 1 
                        if node == t:
                            f = F[P[t] ^ 1]
                            while node != s:
                                eind = P[node]
                                node = V[eind]
                                f = min(f, F[eind ^ 1])
                            node = t
                            while node != s:
                                eind = P[node]
                                node = V[eind]
                                F[eind] += f
                                F[eind ^ 1] -= f
                            total_flow += f
                    else:
                        edge_counter[node] -= 1
                if node == s:
                    break
                node = V[P[node]]
                edge_counter[node] -= 1

        return total_flow

    def calc_flow():
        # Return the (signed) flow in the same order as the edges were created in
        return [(V[eind ^ 1], V[eind], C[eind] - F[eind]) for eind in range(0, len(F), 2)]
    
    return flow_graph(add_edge, Dinic, calc_flow)

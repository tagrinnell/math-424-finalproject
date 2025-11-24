# math-424-finalproject


```
Boruvka's Algorithm (Input: Graph G(V, E))
Output: G' := Minimum Spanning Tree of Input Graph G

    Initialize output graph to have V components

    completed = false
    while not completed
        cheapest_edges = [array of edges]
        for each vertex v in each component:
            for each edge, e of v:
                u_rep = find_set_rep(e.u)
                v_rep = find_set_rep(e.v)

                if u_rep != v_rep and
                        e is the cheapest edge found for component u_rep:
                    cheapest_edges.append(e)


        if size of cheapest_edges > 0:
            for each edge in cheapest_edges:
                Add cheapest_edge to graph G'
                union_set_by_rank(cheapest_edge.u, cheapest_edge.v)
        else :
            completed = true

    return G'
```
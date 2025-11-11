/**
 * Graph implementation using disjoint set (Union find) data structure
 *
 * Sources:
 * https://en.wikipedia.org/wiki/Disjoint-set_data_structure
 * https://www.geeksforgeeks.org/dsa/boruvkas-algorithm-greedy-algo-9/
 * https://www.geeksforgeeks.org/dsa/union-by-rank-and-path-compression-in-union-find-algorithm/
 * https://www.geeksforgeeks.org/dsa/introduction-to-disjoint-set-data-structure-or-union-find-algorithm/
 *
**/
#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "node.hpp"

#include <vector>

class graph {
public:
    int num_vertices;
    std::vector<std::vector<int>> adj_list;

    graph() : num_vertices(0) {}
    graph(int num_vertices) : num_vertices(num_vertices) {}

    void add_edge(int v1, int v2, int weight) {
        adj_list.push_back({v1, v2, weight});
    }

    int find_set_rep(int x) {
        if
    }

private:
};

#endif
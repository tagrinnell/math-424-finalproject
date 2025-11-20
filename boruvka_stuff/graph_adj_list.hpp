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

#include <vector>
#include <array>
#include <unordered_map>
#include <iostream>

class graph_adj_list {
public:
    int num_vertices;
    std::vector<std::vector<int>> edge_list;
    std::vector<int> parents;
    std::vector<int> ranks;

    graph_adj_list() : num_vertices(0) {
        init();
    }


    graph_adj_list(int num_vertices) : num_vertices(num_vertices), parents(num_vertices, -1), ranks(num_vertices) {
        for (int i = 0; i < num_vertices; i++) {
            edge_list.push_back(std::vector<int>(num_vertices, -1));
        }
        init();
    }

    void add_edge(int v1, int v2, int weight) {
        edge_list[v1][v2] = weight;
        edge_list[v2][v1] = weight;
    }

    void add_edge(std::tuple<int, int, int> edge) {
        add_edge(std::get<0>(edge), std::get<1>(edge), std::get<2>(edge));
    }


    int find_set_rep(int x) {
        if (parents[x] == x) {
            return x;
        }
        return find_set_rep(parents[x]);
    }

    void union_set(int x, int y) {
        auto x_rep = find_set_rep(x);
        auto y_rep = find_set_rep(y);

        if (x_rep == y_rep) {
            return;
        }

        if (ranks[x_rep] < ranks[y_rep]) {
            auto tmp = x_rep;
            x_rep = y_rep;
            y_rep = tmp;
        }

        parents[y_rep] = x_rep;
        if (ranks[x_rep] == ranks[y_rep]) {
            ranks[x_rep]++;
        }
    }

    void init(){
        for (int i = 0; i < num_vertices; i++) {
            parents[i] = i;
            ranks[i] = 0;
        }
    }

    void to_string() {
        std::cout << "\n-----------------------------------------------" << std::endl;
        std::cout << "\nGraph to string: " << std::endl;
        for (int i = 0; i < edge_list.size() - 1; i++) {
            for (int j = i + 1; j < edge_list[i].size(); j++) {
                if (edge_list[i][j] != -1) {
                    std::cout << "Weight (" << i << ", " << j << ") = " << edge_list[i][j] << std::endl;
                }
            }
        }
        std::cout << "\n-----------------------------------------------" << std::endl;
    }

    // bool edge_already_exists(int u, int v, int w) {
    //     if (edge_list.find(u) == edge_list.end()) {
    //         return false;
    //     }

    //     for (int i = 0; i < edge_list[u].size(); i++) {
    //         if (std::get<0>(edge_list[u][i]) == v && std::get<1>(edge_list[u][i]) == w) {
    //             return true;
    //         }
    //     }

    //    return false;
    // }
};

#endif

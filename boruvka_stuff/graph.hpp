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
#include <unordered_map>
#include <iostream>

class graph {
public:
    int num_vertices;
    std::vector<std::tuple<int, int, int>> edge_list;
    std::vector<int> parents;
    std::vector<int> ranks;

    graph() : num_vertices(0) {
        init();
    }


    graph(int num_vertices) : num_vertices(num_vertices), parents(num_vertices, -1), ranks(num_vertices) {
        init();
    }

    void add_edge(int v1, int v2, int weight) {
        edge_list.push_back(std::make_tuple(v1, v2, weight));
    }

    /*
    function Find(x) is
        if x.parent != x then
            x.parent := Find(x.parent)
            return x.parent
        else
            return x
            end if
    end function
    */
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
        std::cout << "Graph to string: " << std::endl;
        for (int i = 0; i < edge_list.size(); i++) {
            auto curr_tuple = edge_list[i];
            std::cout << "\t(" << std::get<0>(curr_tuple) << ", " << std::get<1>(curr_tuple) << "), weight = " << std::get<2>(curr_edge) << std::endl;
        }
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

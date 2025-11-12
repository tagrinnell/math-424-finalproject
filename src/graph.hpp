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
    std::unordered_map<int, std::vector<std::tuple<int, int>>> edge_list;
    std::vector<int> parents;
    std::vector<int> ranks;

    graph() : num_vertices(0) {
        init();
    }


    graph(int num_vertices) : num_vertices(num_vertices), parents(num_vertices, -1), ranks(num_vertices) {
        init();
    }

    void add_edge(int v1, int v2, int weight) {
        edge_list[v1].push_back(std::make_tuple(v2, weight));
        edge_list[v2].push_back(std::make_tuple(v1, weight));
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
        // parents(num_vertices);
        // ranks(num_vertices);
        for (int i = 0; i < num_vertices; i++) {
            parents[i] = i;
            ranks[i] = 0;
        }
    }

    void to_string() {
        std::cout << "Graph to string: " << std::endl;
        for (int i = 0; i < num_vertices; i++) {
            auto curr_edge_list = edge_list[i];
            for (int j = 0; j < curr_edge_list.size(); j++) {
                auto curr_edge = curr_edge_list[j];
                std::cout << "\t(" << i << ", " << std::get<0>(curr_edge) << "), weight = " << std::get<1>(curr_edge) << std::endl;
            }
        }
    }

};

#endif

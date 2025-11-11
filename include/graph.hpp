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
#include <unordered_map>

class graph {
public:
    int num_vertices;
    std::unordered_map<int, std::vector<std::tuple<int, int>>> edge_list;
    std::vector<std::vector<node>> tree;


    graph() : num_vertices(0) {}
    graph(int num_vertices) : num_vertices(num_vertices) {
    }

    void add_edge(int v1, int v2, int weight) {
        edge_list[v1].push_back(std::make_tuple(v2, weight));
        edge_list[v2].push_back(std::make_tuple(v1, weight));
    }

    /*
    function Find(x) is
        if x.parent ≠ x then
            x.parent := Find(x.parent)
            return x.parent
        else
            return x
            end if
    end function
    */
    static node* find_set_rep(node* x) {
        if (x->parent != x) {
            return x->parent;
        }
        return find_set_rep(x->parent);
    }

    void union_set(node* x, node* y) {
        node* x_rep = find_set_rep(x);
        node* y_rep = find_set_rep(y);

        if (x_rep == y_rep) {
            return;
        }

        if (x_rep->rank < y_rep->rank) {
            auto tmp = y_rep;
            x_rep = y_rep;
            y_rep = tmp;
        }

        y_rep->parent = x_rep;
        if (x_rep->rank == y_rep->rank) {
            x_rep->rank++;
        }
    }

    void combine_edge_lists(std::vector<std::tuple<int, int>> el_x, std::vector<std::tuple<int, int>> el_y) {

    }

};

#endif

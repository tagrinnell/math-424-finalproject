
#ifndef BORUVKA_HPP
#define BORUVKA_HPP

#include "graph.hpp"

#include <cmath>

// Pseudocode from Wikipedia
/*
algorithm Boruvka is
    input: A weighted undirected graph G = (V, E).
    output: F, a minimum spanning forest of G.

    Initialize a forest F to (V, E′) where E′ = {}.

    completed := false
    while not completed do
        Find the connected components of F and assign to each vertex its component
        Initialize the cheapest edge for each component to "None"
        for each edge uv in E, where u and v are in different components of F:
            let wx be the cheapest edge for the component of u
            if is-preferred-over(uv, wx) then
                Set uv as the cheapest edge for the component of u
            let yz be the cheapest edge for the component of v
            if is-preferred-over(uv, yz) then
                Set uv as the cheapest edge for the component of v
        if all components have cheapest edge set to "None" then
            // no more trees can be merged -- we are finished
            completed := true
        else
            completed := false
            for each component whose cheapest edge is not "None" do
                Add its cheapest edge to E'

function is-preferred-over(edge1, edge2) is
    return (edge2 is "None") or
           (weight(edge1) < weight(edge2)) or
           (weight(edge1) = weight(edge2) and tie-breaking-rule(edge1, edge2))

function tie-breaking-rule(edge1, edge2) is
    The tie-breaking rule; returns true if and only if edge1
    is preferred over edge2 in the case of a tie.

*/
graph boruvka_mst(graph* input_graph) {
    // Construct output graph with N components where N is the number of vertices in the input graph.
    graph output_graph(input_graph->num_vertices);

    bool completed = false;
    while (!completed) {
        std::unordered_map<int, std::tuple<int, int, int>> cheapest_node;

        // Find the cheapest edge for
        for (int i = 0; i < input_graph->edge_list.size(); i++) {
            auto curr_node_edge_list = input_graph->edge_list[i];
            auto x_rep = output_graph.find_set_rep(i);
            for (int j = 0; j < curr_node_edge_list.size(); j++) {
                auto y = std::get<0>(curr_node_edge_list[i]);
                auto y_rep = output_graph.find_set_rep(y);

                // x
                if (x_rep == y_rep) {
                    continue;
                } else {
                    if (cheapest_node.find(x_rep) == cheapest_node.end() ||
                            std::get<2>(cheapest_node[i]) > std::get<1>(curr_node_edge_list[i])) {
                        cheapest_node[i] = std::make_tuple(i, y, std::get<1>(curr_node_edge_list[i]));
                    }

                }

            }
        }

    }

    return output_graph;
}

#endif

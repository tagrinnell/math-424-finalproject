#ifndef BORUVKA_HPP
#define BORUVKA_HPP

#include "graph.hpp"

#include <iostream>
#include <cmath>
#ifdef _OPENMP
    #include <omp.h>
#endif

inline bool tie_breaking_rule(std::tuple<int, int, int> edge1, std::tuple<int, int, int> edge2);
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
graph boruvka_mst(graph input_graph) {
    graph output_graph(input_graph.num_vertices);

    std::unordered_map<int, std::vector<int>> component_list;
    std::unordered_map<int, std::tuple<int, int, int>> cheapest_node;
    auto input_edge_list = input_graph.edge_list;

    bool completed = false;
    while (!completed) {
        std::cout << "New Iteration" << std::endl;
        // Find all of the components
        for (int i = 0; i < output_graph.num_vertices; i++) {
            auto x = output_graph.find_set_rep(i);
            // This component hasn't been found yet, create a new component in the map
            if (component_list.find(x) == component_list.end()) {
                component_list[x] = {i};
            } else {
                component_list[x].push_back(i); // Add to the component list
            }
        }

        // Initialize the cheapest edge for each component to None
        cheapest_node.clear();

        for (int i = 0; i < input_edge_list.size(); i++) {
            auto u = std::get<0>(input_edge_list[i]);
            auto v = std::get<1>(input_edge_list[i]);
            auto w = std::get<2>(input_edge_list[i]);

            auto u_rep = output_graph.find_set_rep(u);
            auto v_rep = output_graph.find_set_rep(v);

            // u and v belong to different components
            if (u_rep != v_rep) {

                // u's component doesn't currently have a cheapest node or
                // the weight of the current edge is smaller than the current smallest edge weight
                if (cheapest_node.find(u_rep) == cheapest_node.end() ||
                      w < std::get<2>(cheapest_node[u_rep]) ||
                      tie_breaking_rule(input_edge_list[i], cheapest_node[u_rep])) {
                    cheapest_node[u_rep] = input_edge_list[i];
                }

                if (cheapest_node.find(v_rep) == cheapest_node.end() ||
                      w < std::get<2>(cheapest_node[v_rep]) ||
                      tie_breaking_rule(input_edge_list[i], cheapest_node[v_rep])) {
                    cheapest_node[v_rep] = input_edge_list[i];
                }

            }
            // else {
                // As an optimization, one could remove from G each edge that is found to connect two
                // vertices in the same component, so that it does not contribute to the time for searching
                // for cheapest edges in later components.
           // }

        }

        if (cheapest_node.size() != 0) {
            for (const auto& x : cheapest_node) {
                std::cout << "  Adding edge from " << std::get<0>(x.second) << " to " << std::get<1>(x.second) << " with weight " << std::get<2>(x.second) << std::endl;
                output_graph.add_edge(x.second);
                output_graph.union_set(std::get<0>(x.second), std::get<1>(x.second));
            }
        } else {
            completed = true;
        }

        component_list.clear();
    }

    return output_graph;
}

// The tie-breaking rule; returns true if and only if edge1
// is preferred over edge2 in the case of a tie.
// A tie-breaking rule which orders edges first by source, then by destination,
// will prevent creation of a cycle, resulting in the minimal spanning tree {ab, bc}.
inline bool tie_breaking_rule(std::tuple<int, int, int> edge1, std::tuple<int, int, int> edge2) {
    return (std::get<2>(edge1) == std::get<2>(edge2)) && (std::get<0>(edge1) < std::get<0>(edge2)) && (std::get<1>(edge1) < std::get<1>(edge2));
}

/*
graph boruvka_mst_old(graph* input_graph) {
    // Construct output graph with N components where N is the number of vertices in the input graph.
    graph output_graph(input_graph->num_vertices);

    bool completed = false;
    std::unordered_map<int, std::tuple<int, int, int>> cheapest_node;

    // Initialize component list to contain each of the vertices
    std::vector<std::vector<int>> component_list;
    for (int i = 0; i < input_graph->num_vertices; i++) {
        component_list.push_back({i});
    }
    while (!completed) {
        std::cout << "New iteration" << std::endl;
        cheapest_node.clear();

        // For each component
        for (int i = 0; i < component_list.size(); i++) {
            auto cur

        }
        // Find the cheapest edge for
        for (int i = 0; i < input_graph->edge_list.size(); i++) {
            auto curr_node_edge_list = input_graph->edge_list[i];
            auto x_rep = output_graph.find_set_rep(i);
            for (int j = 0; j < curr_node_edge_list.size(); j++) {
                auto y = std::get<0>(curr_node_edge_list[j]);
                auto y_rep = output_graph.find_set_rep(y);

                // x and y belong to the same set
                if (x_rep == y_rep) {
                    continue;
                } else {            // x_rep and y_rep
                    if ((cheapest_node.find(x_rep) == cheapest_node.end()) ||                                  // The cheapest node hasn't been found yet
                            // cheapest node has a 3-tuple //node_edge list is a map of 2-tuples
                            (std::get<2>(cheapest_node[x_rep]) > std::get<1>(curr_node_edge_list[j]))         // Current node has a smaller weight than the current cheapest edge
                        ) {
                        cheapest_node[x_rep] = std::make_tuple(i, y, std::get<1>(curr_node_edge_list[j]));
                    }

                }
            }

        }

        auto num_cheapest_found = 0;
        for (int i = 0; i < output_graph.num_vertices; i++) {
            // We have a cheapest node.  Add to output_graph
            // increment count of found cheapest edges.
            if (cheapest_node.find(i) != cheapest_node.end()// &&
                // (output_graph.edge_list.find(i) == output_graph.edge_list.end())            //    std::get<1>(output_graph.edge_list.find(i))
                ) {
                auto u = std::get<0>(cheapest_node[i]);
                auto v = std::get<1>(cheapest_node[i]);
                auto w = std::get<2>(cheapest_node[i]);

                if (!output_graph.edge_already_exists(u, v, w)) {
                        std::cout << "\tAdding edge from " << u << " to " << v << " weighing " << w << std::endl;
                        output_graph.add_edge(u, v, w);
                        output_graph.union_set(i, v);   // Join vertex v to i (representative of this tree)
                        num_cheapest_found++;
               }
            }
        }
        if (num_cheapest_found == 0) {
            completed = true;
        }
    }

    return output_graph;
}*/

#endif
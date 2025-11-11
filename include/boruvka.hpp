
#ifndef BORUVKA_HPP
#define BORUVKA_HPP

#include "graph.hpp"
#include "node.hpp"

#include <cmath>

graph boruvka_mst(graph* input_graph) {
    graph output_graph(input_graph->num_vertices);

    // Construct output graph with N components where N is the number of vertices in the input graph.
    for (int i = 0; i < output_graph.num_vertices; i++) {
        node new_node;
        new_node.parent = &new_node;
        new_node.node_num = i;
        new_node.rank = 0;
        new_node.cheapest_edge = std::make_tuple(i, i, RAND_MAX);

        output_graph.vertices.push_back({new_node});
    }


    bool completed = false;
    while (!completed) {
        std::vector<int> cheapest_node = {0, 0, RAND_MAX};


    }

    return output_graph;
}

#endif

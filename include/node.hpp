/**
 * Backing node structure for the graph
 */

#include <tuple>
#include <cmath>

class node {
public:
    node* parent;
    node* next_node;
    int node_num;
    int rank;
    std::tuple<int, int, int> cheapest_edge;

    node() : node_num(0) {}
    node(node* parent, node* next_node, int node_num) : node_num(node_num), parent(parent), next_node(next_node), rank(0) {
        cheapest_edge = std::make_tuple(node_num, node_num, RAND_MAX);
    }
};

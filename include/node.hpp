/**
 * Backing node structure for the graph
 */

class node {
public:
    node* parent;
    node* next_node;
    int node_num;

    node() : node_num(0) {}
    node(node* parent, node* next_node, int node_num) : node_num(node_num), parent(parent), next_node(next_node) {}

};
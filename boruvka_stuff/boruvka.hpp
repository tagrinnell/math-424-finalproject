#ifndef BORUVKA_HPP
#define BORUVKA_HPP

#include "graph_adj_list.hpp"

#include <iostream>
#include <cmath>
#include <chrono>

// #define _MPI true

#ifdef _OPENMP
    #include <omp.h>
#endif
// #ifdef _MPI
    #include <mpi.h>
// #endif

#define ROOT_RANK 0

#define PARENT_VEC_TAG 0
#define SCATTER_TAG 1
#define REMAINDER_TAG 2
#define EDGE_TAG 3

inline bool tie_breaking_rule(int u, int v, int w, std::tuple<int, int, int> cheapest_node_edge);
inline bool tie_breaking_rule(int u, int v, int w, int i, int j, int k);
inline void print_tuple(std::tuple<int, int, int, bool> edge);

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
graph_adj_list boruvka_mst(graph_adj_list input_graph) {
    graph_adj_list output_graph(input_graph.num_vertices);

    // Component representative -> Edge(i, j, w)
    std::unordered_map<int, std::tuple<int, int, int>> cheapest_node;
    auto input_edge_list = input_graph.edge_list;

    bool completed = false;
    while (!completed) {
        std::cout << "New Iteration" << std::endl;
        // Initialize the cheapest edge for each component to None
        cheapest_node.clear();

        for (int i = 0; i < input_edge_list.size(); i++) {
            for (int j = 0; j < input_edge_list[i].size(); j++) {
                if (i != j) {
                    auto u = i;
                    auto v = j;
                    auto w = input_edge_list[i][j];

                    if (w > 0) {
                        auto u_rep = output_graph.find_set_rep(u);
                        auto v_rep = output_graph.find_set_rep(v);

                        // u and v belong to different components
                        if (u_rep != v_rep) {

                            // u's component doesn't currently have a cheapest node or
                            // the weight of the current edge is smaller than the current smallest edge weight
                            if (cheapest_node.find(u_rep) == cheapest_node.end() ||
                                    w < std::get<2>(cheapest_node[u_rep]) ||
                                    tie_breaking_rule(u, v, w, cheapest_node[u_rep])) {
                                cheapest_node[u_rep] = std::make_tuple(u, v, w);
                            }

                            // if (cheapest_node.find(v_rep) == cheapest_node.end() ||
                            //         w < std::get<2>(cheapest_node[v_rep]) ||
                            //         tie_breaking_rule(u, v, w, cheapest_node[v_rep])) {
                            //     cheapest_node[v_rep] = std::make_tuple(u, v, w);
                            // }

                        }
                    }
                }
            }
        }

        if (cheapest_node.size() != 0) {
            for (auto& x : cheapest_node) {
                std::cout << "  Adding edge from " << std::get<0>(x.second) << " to " << std::get<1>(x.second) << " with weight " << int(std::get<2>(x.second)) << std::endl;
                output_graph.add_edge(x.second);
                output_graph.union_set(std::get<0>(x.second), std::get<1>(x.second));
            }
        } else {
            completed = true;
        }

    }

    return output_graph;
}

// is preferred over edge2 in the case of a tie.
// A tie-breaking rule which orders edges first by source, then by destination,
// will prevent creation of a cycle, resulting in the minimal spanning tree {ab, bc}.
inline bool tie_breaking_rule(int u, int v, int w, std::tuple<int, int, int> cheapest_node_edge) {
    return (w == std::get<2>(cheapest_node_edge)) && (u < std::get<0>(cheapest_node_edge)) && (v < std::get<1>(cheapest_node_edge));
}

inline void print_tuple(std::tuple<int, int, int, bool> edge) {
    std::cout << "\tEdge from " << std::get<0>(edge) << " to " << std::get<1>(edge) << " with weight " << std::get<2>(edge) << ", visited = " << std::get<3>(edge) << std::endl;
}

#ifdef _OMP
graph boruvka_mst_openmp(graph input_graph, int num_threads) {
    graph_adj_list output_graph(input_graph.num_vertices);

    // Component representative -> Edge(i, j, w)
    std::unordered_map<int, std::tuple<int, int, int>> cheapest_node;
    auto input_edge_list = input_graph.edge_list;

    bool completed = false;
    while (!completed) {
        std::cout << "New Iteration" << std::endl;
        // Initialize the cheapest edge for each component to None
        cheapest_node.clear();

        #pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < input_edge_list.size(); i++) {
            for (int j = 0; j < input_edge_list[i].size(); j++) {
                if (i != j) {
                    auto u = i;
                    auto v = j;
                    auto w = input_edge_list[i][j];

                    if (w > 0) {
                        auto u_rep = output_graph.find_set_rep(u);
                        auto v_rep = output_graph.find_set_rep(v);

                        // u and v belong to different components
                        if (u_rep != v_rep) {
                            // u's component doesn't currently have a cheapest node or
                            // the weight of the current edge is smaller than the current smallest edge weight
                            if (cheapest_node.find(u_rep) == cheapest_node.end() ||
                                    w < std::get<2>(cheapest_node[u_rep]) ||
                                    tie_breaking_rule(u, v, w, cheapest_node[u_rep])) {
                                cheapest_node[u_rep] = std::make_tuple(u, v, w);
                            }
                        }
                    }
                }
            }
        }

        if (cheapest_node.size() != 0) {
            for (auto& x : cheapest_node) {
                std::cout << "  Adding edge from " << std::get<0>(x.second) << " to " << std::get<1>(x.second) << " with weight " << int(std::get<2>(x.second)) << std::endl;
                output_graph.add_edge(x.second);
                output_graph.union_set(std::get<0>(x.second), std::get<1>(x.second));
            }
        } else {
            completed = true;
        }

    }

    return output_graph;
}
#endif

#ifndef _OPENMP && _MPI

// TODO change to receive
std::vector<int> parent_arr_receive(int size, int rank, int parent_vector_size) {
    std::vector<int> r_vector(parent_vector_size);
    // Receive the Canonical parent vector from root node
    if (rank != 0) {
        MPI_Bcast(&r_vector[0], parent_vector_size, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    }
    return r_vector;
}

void parent_arr_send(int size, int rank, std::vector<int> parent_vector) {
    if (rank == 0) {
        MPI_Bcast(&parent_vector[0], parent_vector.size(), MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    } else {
        std::cout << "Sending from rank " << rank << std::endl;
        MPI_Send(&parent_vector[0], parent_vector.size(), MPI_INT, rank, PARENT_VEC_TAG, MPI_COMM_WORLD);
    }
}

std::vector<int> edge_receive(int size, int rank, int dispatched_components) {
    std::vector<int> receive_edges;
    // Receiving edges from
    if (rank == 0) {
        for (int i = 1; i < size || i < dispatched_components; i++) {
            // Receive number of edges going to be transmitted by non-root rank processes
            int num_rec;
            MPI_Recv(&num_rec, 1, MPI_INT, i, EDGE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int old_size = receive_edges.size();
            receive_edges.resize(old_size + num_rec);
            MPI_Recv(&receive_edges[old_size], num_rec, MPI_INT, i, EDGE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    return receive_edges;
}

void edge_send(int size, int rank, std::vector<int> edges_to_send) {
    // Send edges_to_send to the root node
    if (rank != 0) {
        MPI_Send(&edges_to_send[0], edges_to_send.size(), MPI_INT, ROOT_RANK, EDGE_TAG, MPI_COMM_WORLD);
    }
}


// TODO test
std::unordered_map<int, bool> component_arr_receive(int size, int rank) {
    std::unordered_map<int, bool> ret_map;

    // Rank 0 --> receiving from all other processes
    if (rank != 0) {
        // Receive # of components to receive
        int receive_size;
        MPI_Bcast(&receive_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::vector<int> receiving_vector(receive_size);
        MPI_Scatter(&receiving_vector[0], 0, MPI_INT, &receiving_vector[0], receive_size, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

        if (rank == size - 1) {
            int scatter_remainder;
            MPI_Recv(&scatter_remainder, 1, MPI_INT, ROOT_RANK, SCATTER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (scatter_remainder > 0) {
                receiving_vector.resize(receive_size + scatter_remainder);
                MPI_Recv(&receiving_vector[receiving_vector.size() - scatter_remainder], scatter_remainder, MPI_INT, ROOT_RANK, REMAINDER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        // initialize map to return for easy searching
        for (auto i : receiving_vector) {
            ret_map[i] = true;
        }
    }
    return ret_map;
}

// TODO Test
// TODO Adjust for case where component_list.size() < # available processes = size - 1
void component_arr_send(int size, int rank, std::vector<int> component_list) {
    if (rank == 0) {
        // We can scatter multiple components per process
        if (component_list.size() > size - 1) {
            // Scatter count
            auto scatter_count = component_list.size() / (size- 1);  // Divide among the ranks other than 1
            auto scatter_remainder = component_list.size() - (size - 1) * scatter_count;

            MPI_Bcast(&scatter_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // send_count --> how much is sent to each process
            MPI_Scatter(&component_list[0], scatter_count, MPI_INT, NULL, 0, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

            // Send the number of remaining elements to the last process
            MPI_Send(&scatter_remainder, 1, MPI_INT, size - 1, REMAINDER_TAG, MPI_COMM_WORLD);
            MPI_Send(&component_list[component_list.size() - scatter_remainder], scatter_remainder, MPI_INT, size- 1, REMAINDER_TAG, MPI_COMM_WORLD);
        } else {
            // TODO Send
        }
    } else {
        std::cout << "Sending from rank " << rank << std::endl;
        // MPI_Send(&parent_vector[0], parent_vector.size(), MPI_INT, 0, MPI_PARENT_VEC_TAG, MPI_COMM_WORLD);
    }
}

inline bool tie_breaking_rule(int u, int v, int w, int i, int j, int k) {
    return (w == k) && (u < i) && (v < j);
}

graph_adj_list boruvka_mst_mpi(graph_adj_list input_graph) {
    // Idea:
    // All processes will have the list of edges

    /*
    while (not completed)
        Broadcast: Parent vector (Or changes to the parent vector)
        Scatter: Component List (Representative)

        for each component in the component list
            find cheapest edge for each component (but only that component)
            union set and

        Gather:
            Added edges
            Parent vectors (Need to )
    */
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size <= 2) {
        return boruvka_mst(input_graph);
    }

    auto n_vertices = input_graph.num_vertices;
    auto input_edge_list = input_graph.edge_list;
    graph_adj_list output_graph(n_vertices);

    bool completed = false;
    while (!completed) {
        int dispatched_components = 0;
        // Find Components
        if (rank == 0) {
            std::vector<int> component_list;
            std::unordered_map<int, bool> comp_test;
            for (int i = 0; i < n_vertices; i++) {
                auto rep = output_graph.find_set_rep(i);
                if (comp_test.find(rep) == comp_test.end()) {
                    comp_test[rep] = true;
                    component_list.push_back(rep);
                }
            }
            if (component_list.size() < size - 1) {
                dispatched_components = component_list.size();
            }
            // Scatter components across the world
            component_arr_send(size, rank, component_list);
            parent_arr_send(size, rank, output_graph.parents);
        } else {
            // Receive component number
            auto comp_map = component_arr_receive(size, rank);
            if (comp_map.size() == 0) {
                completed = true;
                break;
            }
            output_graph.parents = parent_arr_receive(size, rank, output_graph.num_vertices);
            std::unordered_map<int, std::tuple<int, int, int>> cheapest_node;

            for (int i = 0; i < n_vertices; i++) {
                auto find_rep = output_graph.find_set_rep(i);

                // We have a component that was sent to this process
                if (comp_map.find(find_rep) != comp_map.end()) {
                    for (int j = 0; j < input_edge_list[i].size(); j++) {
                        auto u = i;
                        auto v = j;
                        auto w = input_edge_list[i][j];

                        if (w > 0) {
                            auto u_rep = find_rep;
                            auto v_rep = output_graph.find_set_rep(v);

                            // u and v belong to different components
                            if (u_rep != v_rep) {

                                // Same algorithm except we only check the
                                if (cheapest_node.find(u_rep) == cheapest_node.end() ||
                                        w < std::get<2>(cheapest_node[u_rep]) ||
                                        tie_breaking_rule(u, v, w, cheapest_node[u_rep])) {
                                    cheapest_node[u_rep] = std::make_tuple(u, v, w);
                                }

                                // if (cheapest_node.find(v_rep) == cheapest_node.end() ||
                                //         w < std::get<2>(cheapest_node[v_rep]) ||
                                //         tie_breaking_rule(u, v, w, cheapest_node[v_rep])) {
                                //     cheapest_node[v_rep] = std::make_tuple(u, v, w);
                                // }

                            }
                        }
                    }
                }
            }

            // TODO
            // Transmit all found edges back to the rank = 0
            std::vector<int> send_vec;
            for (auto x : cheapest_node) {
                send_vec.push_back(std::get<0>(x.second));
                send_vec.push_back(std::get<1>(x.second));
                send_vec.push_back(std::get<2>(x.second));
            }
            edge_send(size, rank, send_vec);
        }

        // Find cheapest edge for all components
        // Using each node that is a part of the
        if (rank == 0) {
            // Gather all of the edges
            // Add these to the mst and perform set unions
            auto received_edges = edge_receive(size, rank, dispatched_components);
            // Split into 3-strides --> 0, 1, 2 -> edge(0, 1) = 2
            for (int i = 0; i < received_edges.size() - 2; i = i + 3) {
                output_graph.add_edge(received_edges[i], received_edges[i + 1], received_edges[i + 2]);
                output_graph.union_set(received_edges[i], received_edges[i + 1]);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    return output_graph;
}

// https://github.com/nikitawani07/MST-Parallel/blob/master/src/boruvka.c
graph_adj_list mpi_wrapper(graph_adj_list input_graph, int  argc, char** argv) {
    auto err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        std::cout << "MPI failed somehow" << std::endl;
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
    auto out = boruvka_mst_mpi(input_graph);
    const auto t2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << "Serial Implementation ran for " << ms.count() << std::endl;

    MPI_Finalize();

    return out;
}

// #endif

#endif

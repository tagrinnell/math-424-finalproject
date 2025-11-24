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
#define COMPLETED_TAG 4

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

#ifndef _OPENMP

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
        // std::cout << "Root entering parent_arr_send" << std::endl;
        MPI_Bcast(&parent_vector[0], parent_vector.size(), MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    } else {
        // std::cout << "Sending from rank " << rank << std::endl;
        MPI_Send(&parent_vector[0], parent_vector.size(), MPI_INT, rank, PARENT_VEC_TAG, MPI_COMM_WORLD);
    }
}

std::vector<int> edge_receive(int size, int rank, int dispatched_components) {
    std::vector<int> receive_edges;
    // Receiving edges from
    if (rank == 0) {
        // std::cout << "Root entering edge_receive" << std::endl;
        for (int i = 1; i < size && i < dispatched_components; i++) {
            // Receive number of edges going to be transmitted by non-root rank processes
            int num_rec = 0;
            MPI_Recv(&num_rec, 1, MPI_INT, i, EDGE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int old_size = receive_edges.size();
            receive_edges.resize(old_size + num_rec);
            // std::cout << "Old size " << old_size << "; New size " << receive_edges.size() << std::endl;
            MPI_Recv(&receive_edges[old_size], num_rec, MPI_INT, i, EDGE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // for (auto i : receive_edges) {
            //     std::cout << i << " ";
            // }
            // std::cout << std::endl;
        }
    }
    return receive_edges;
}

void edge_send(int size, int rank, std::vector<int> edges_to_send) {
    // Send edges_to_send to the root node
    if (rank != 0) {
        int size_send = edges_to_send.size();
        MPI_Send(&size_send, 1, MPI_INT, ROOT_RANK, EDGE_TAG, MPI_COMM_WORLD);
        MPI_Send(&edges_to_send[0], edges_to_send.size(), MPI_INT, ROOT_RANK, EDGE_TAG, MPI_COMM_WORLD);
    }
}

// TODO test
std::unordered_map<int, bool> component_arr_scatter(int size, int rank, std::vector<int> component_list = std::vector<int>()) {
    std::unordered_map<int, bool> ret_map;
    std::vector<int> receiving_vector;

    // Rank 0 --> receiving from all other processes
    if (rank == 0) {
        // std::cout << "Root entering component_arr_scatter" << std::endl;
         // We can scatter multiple components per process
        if (component_list.size() > size) {
            // std::cout << " Input scatter map: ";
            // for (auto y : component_list) {
            //     std::cout << y << " ";
            // }
            // std::cout << std::endl;
            // Scatter count
            auto scatter_count = component_list.size() / size;  // Divide among the ranks other than 1
            auto scatter_remainder = component_list.size() - size * scatter_count;

            // std::cout << "Broadcasting" << std::endl;;
            MPI_Bcast(&scatter_count, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

            // std::cout << "Scattering" << std::endl;;
            receiving_vector.resize(scatter_count);
            // send_count --> how much is sent to each process
            MPI_Scatter(&component_list[0], scatter_count, MPI_INT,
                        &receiving_vector[0], scatter_count, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

            // std::cout << "Scattering remainder " << scatter_remainder << std::endl;

            // Send the number of remaining elements to the last process
            MPI_Send(&scatter_remainder, 1, MPI_INT, size - 1, REMAINDER_TAG, MPI_COMM_WORLD);
            if (scatter_remainder > 0) {
                MPI_Send(&component_list[component_list.size() - scatter_remainder], scatter_remainder, MPI_INT, size- 1, REMAINDER_TAG, MPI_COMM_WORLD);
            }
        } else if (component_list.size() > 1) {
            // TODO Scatter a single component to each of the remaining processes
            int zero = 0;
            int terminate_char = -1;
            MPI_Bcast(&zero, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
            for (int i = 1; i < size; i++) {
                if (i < component_list.size()) {
                    MPI_Send(&component_list[i], 1, MPI_INT, i, SCATTER_TAG, MPI_COMM_WORLD);
                } else {
                    MPI_Send(&terminate_char, 1, MPI_INT, i, SCATTER_TAG, MPI_COMM_WORLD);
                }
            }

            // Send termination value to the rest of the processes
            for (int i = component_list.size(); i < size; i++) {
            }
            receiving_vector.push_back(component_list[0]);
        } else {

        }
    } else {
        // std::cout << "Rank " << rank << " receiving components" << std::endl;
        // Receive # of components to receive
        int receive_size;
        MPI_Bcast(&receive_size, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

        // Special Code to exit
        if (receive_size == -1) {
            receiving_vector.push_back(-1);
        } else if (receive_size == 0) {
            receiving_vector.resize(1);
            MPI_Recv(&receiving_vector[0], 1, MPI_INT, ROOT_RANK, SCATTER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            receiving_vector.resize(receive_size);
            std::vector<int> component_list(receive_size);
            // std::cout << "  Receive size " << receive_size << std::endl;
            MPI_Scatter(&component_list[0], receive_size, MPI_INT,
                &receiving_vector[0], receive_size, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

            // std::cout << "\n  Rank " << rank << " after scatter " << std::endl;
            // for (auto i : receiving_vector) {
            //     std::cout << i << " ";
            // }

            if (rank == size - 1) {
                // std::cout << "Reaching scatter remainder section " << std::endl;
                int scatter_remainder;
                MPI_Recv(&scatter_remainder, 1, MPI_INT, ROOT_RANK, REMAINDER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // std::cout << "Scatter reminder " << scatter_remainder << std::endl;
                if (scatter_remainder > 0) {
                    receiving_vector.resize(receive_size + scatter_remainder);
                    MPI_Recv(&receiving_vector[receiving_vector.size() - scatter_remainder], scatter_remainder, MPI_INT, ROOT_RANK, REMAINDER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

    }
    // std::cout << "\t Rank " << rank << ": ";
    // initialize map to return for easy searching
    for (auto i : receiving_vector) {
        // std::cout << i << " ";
        ret_map[i] = true;
    }
    // std::cout << std::endl;
    return ret_map;
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

    // std::cout << "Rank: " << rank << " - Size: " << size << std::endl;

    if (size == 1) {
        return boruvka_mst(input_graph);
    }

    auto n_vertices = input_graph.num_vertices;
    auto input_edge_list = input_graph.edge_list;
    graph_adj_list output_graph(n_vertices);

    bool completed = false;
    bool early_exit = false;
    while (!completed) {
        int dispatched_components = 0;
        std::unordered_map<int, bool> comp_map;
        // Find Components
        if (rank == 0) {
            // std::cout << "\nNew Iteration" << std::endl;
            // std::cout << "Completed: " << completed << std::endl;
            std::vector<int> component_list;
            std::unordered_map<int, bool> comp_test;
            for (int i = 0; i < n_vertices; i++) {
                auto rep = output_graph.find_set_rep(i);
                if (comp_test.find(rep) == comp_test.end()) {
                    comp_test[rep] = true;
                    component_list.push_back(rep);
                }
            }
            // std::cout << "\tComponent_list size: " << component_list.size() << std::endl;

            if (component_list.size() < size) {
                dispatched_components = component_list.size();
            } else {
                dispatched_components = size;
            }

            if (component_list.size() == 1) {
                completed = true;
                // Send nonsense to the rest of the nodes that exist
                int termination = -1;
                MPI_Bcast(&termination, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
                break;
            }
            // MPI_Bcast(&completed, 1, MPI_CXX_BOOL, ROOT_RANK, MPI_COMM_WORLD);
            // if (completed) {
            //     break;
            // }
            comp_map = component_arr_scatter(size, rank, component_list);
            parent_arr_send(size, rank, output_graph.parents);
        } else {
            // Check for valid Component map
            // MPI_Bcast(&completed, 1, MPI_CXX_BOOL, ROOT_RANK, MPI_COMM_WORLD);
            comp_map = component_arr_scatter(size, rank);
            if (comp_map.find(-1) != comp_map.end()) {
                // early_exit = true;
                break;
            } else {
                output_graph.parents = parent_arr_receive(size, rank, output_graph.num_vertices);
            }
        }


        // Receive component number
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

        if (rank != 0) {
            // std::cout << "Rank: " << rank<< " ";
            // Transmit all found edges back to the rank = 0
            std::vector<int> send_vec;
            for (auto x : cheapest_node) {
                send_vec.push_back(std::get<0>(x.second));
                send_vec.push_back(std::get<1>(x.second));
                send_vec.push_back(std::get<2>(x.second));
            }

            edge_send(size, rank, send_vec);
            // MPI_Bcast(&completed, 1, MPI_CXX_BOOL, ROOT_RANK, MPI_COMM_WORLD);

        }
        if (rank == 0) {
            for (auto x : cheapest_node) {
                output_graph.add_edge(std::get<0>(x.second), std::get<1>(x.second), std::get<2>(x.second));
                output_graph.union_set(std::get<0>(x.second), std::get<1>(x.second));
            }

            auto received_edges = edge_receive(size, rank, dispatched_components);

            // Split into 3-strides --> 0, 1, 2 -> edge(0, 1) = 2
            for (int i = 0; i < received_edges.size() - 1; i = i + 3) {
                // std::cout << "  Adding edge from " << received_edges[i] << " to " << received_edges[i + 1] << " with weight " << int(received_edges[i + 2]) << std::endl;
                output_graph.add_edge(received_edges[i], received_edges[i + 1], received_edges[i + 2]);
                output_graph.union_set(received_edges[i], received_edges[i + 1]);
            }

        }
    }
    return output_graph;
}

bool assert_correctness(graph_adj_list a, graph_adj_list b) {
    if (a.edge_list.size() != b.edge_list.size() && a.edge_list[0].size() != b.edge_list[0].size()) {
        return false;
    }
    for (int i = 0; i < a.edge_list.size(); i++) {
        for (int j = 0; j < a.edge_list[0].size(); j++) {
            if (a.edge_list[i][j] != b.edge_list[i][j]) {
                return false;
            }
        }
    }
    return true;
}

// https://github.com/nikitawani07/MST-Parallel/blob/master/src/boruvka.c
graph_adj_list mpi_wrapper(graph_adj_list input_graph, int  argc, char** argv) {
    std::chrono::_V2::system_clock::time_point t1, t2;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    graph_adj_list baseline;
    std::chrono::duration<double, std::milli> mpi_ms, serial_ms;
    // std::ofstream output_stats("output_results.log", std::ios::app);
    // std::ofstream output_stats("output_results.log", std::ios::app | std::ios::out | std::ios::in);

    if (rank == 0) {
        std::clog << "\nBaseline " << std::endl;
        // Run Baseline first()
        // std::chrono::_V2::system_clock::time_point t1, t2;
        auto t1 = std::chrono::high_resolution_clock::now();
        baseline = boruvka_mst(input_graph);
        auto t2 = std::chrono::high_resolution_clock::now();

        serial_ms = t2 - t1;

        std::clog << "Serial Implementation ran for " << serial_ms.count() << std::endl;

    }

    MPI_Barrier(MPI_COMM_WORLD);

    t1 = std::chrono::high_resolution_clock::now();
    auto out = boruvka_mst_mpi(input_graph);
    t2 = std::chrono::high_resolution_clock::now();

    if (rank == 0) {
        // std::cout << "Curr Rank " << rank << std::endl;
        // out.to_string();
        std::clog << "Checking for correctness... ";
        if (assert_correctness(baseline, out)) {
            std::clog << "Correct!!" << std::endl;
        } else {
            std::clog << "Something went wrong :(" << std::endl;
        }

        mpi_ms = t2 - t1;
        std::clog  << "MPI Implementation ran for " << mpi_ms.count() << std::endl;
        std::cout << serial_ms.count() << ", " << mpi_ms.count() << std::endl;
    }

    // MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Finalize();

    return out;
}

#endif

#endif

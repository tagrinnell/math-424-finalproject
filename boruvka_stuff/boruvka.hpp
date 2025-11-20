#ifndef BORUVKA_HPP
#define BORUVKA_HPP

#include "graph_adj_list.hpp"

#include <iostream>
#include <cmath>
#include <chrono>

#ifdef _OPENMP
    #include <omp.h>
#endif
// #ifdef _MPI
    #include <mpi.h>
// #endif

#define MPI_PARENT_VEC_TAG 0

inline bool tie_breaking_rule(int u, int v, int w, std::tuple<int, int, int> cheapest_node_edge);
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

                            if (cheapest_node.find(v_rep) == cheapest_node.end() ||
                                    w < std::get<2>(cheapest_node[v_rep]) ||
                                    tie_breaking_rule(u, v, w, cheapest_node[v_rep])) {
                                cheapest_node[v_rep] = std::make_tuple(u, v, w);
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

        // component_list.clear();
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

// graph boruvka_mst_openmp(graph input_graph, int num_threads) {
//     graph output_graph(input_graph.num_vertices);

//     std::unordered_map<int, std::tuple<int, int, int, bool>> cheapest_node;
//     auto input_edge_list = input_graph.edge_list;

//     bool completed = false;
//     while (!completed) {
//         std::cout << "New Iteration" << std::endl;
//         // Initialize the cheapest edge for each component to None
//         cheapest_node.clear();

//         #pragma omp parallel for num_threads(num_threads)
//         for (int i = 0; i < input_edge_list.size(); i++) {

//             auto curr_edge = input_edge_list[i];
//             auto u = std::get<0>(curr_edge);
//             auto v = std::get<1>(curr_edge);
//             auto w = std::get<2>(curr_edge);

//             auto u_rep = output_graph.find_set_rep(u);
//             auto v_rep = output_graph.find_set_rep(v);

//             // u and v belong to different components
//             if (u_rep != v_rep) {

//                 // u's component doesn't currently have a cheapest node or
//                 // the weight of the current edge is smaller than the current smallest edge weight
//                 #pragma omp critical
//                 {
//                     if (cheapest_node.find(u_rep) == cheapest_node.end() ||
//                         w < std::get<2>(cheapest_node[u_rep]) ||
//                         tie_breaking_rule(input_edge_list[i], cheapest_node[u_rep])) {
//                         if (cheapest_node.find(u_rep) != cheapest_node.end()) {
//                             std::get<3>(cheapest_node[u_rep]) = false;
//                         }
//                         cheapest_node[u_rep] = input_edge_list[i];
//                         // std::get<3>(input_edge_list[i]) = true;
//                     }
//                 }

//                 #pragma omp critical
//                 {
//                     if (cheapest_node.find(v_rep) == cheapest_node.end() ||
//                         w < std::get<2>(cheapest_node[v_rep]) ||
//                         tie_breaking_rule(input_edge_list[i], cheapest_node[v_rep])) {
//                         if (cheapest_node.find(v_rep) != cheapest_node.end()) {
//                             std::get<3>(cheapest_node[v_rep]) = false;
//                         }
//                         cheapest_node[v_rep] = input_edge_list[i];
//                         // std::get<3>(input_edge_list[i]) = true;
//                     }
//                 }


//             }
//             //     else {
//             //         // As an optimization, one could remove from G each edge that is found to connect two
//             //         // vertices in the same component, so that it does not contribute to the time for searching
//             //         // for cheapest edges in later components.
//             //         // std::cout << "Removing edge (connects within a single component): ";
//             //         // print_tuple(curr_edge);
//             //         input_edge_list.erase(input_edge_list.begin() + i);
//             //         i--;
//             //    }
//             // }

//         }

//         if (cheapest_node.size() != 0) {
//             for (const auto& x : cheapest_node) {
//                 #pragma omp critical
//                 {
//                     std::cout << "  Adding edge from " << std::get<0>(x.second) << " to " << std::get<1>(x.second) << " with weight " << std::get<2>(x.second) << std::endl;
//                     output_graph.add_edge(x.second);
//                     output_graph.union_set(std::get<0>(x.second), std::get<1>(x.second));
//                 }
//             }
//         } else {
//             completed = true;
//         }

//     }

//     return output_graph;
// }


// https://github.com/nikitawani07/MST-Parallel/blob/master/src/boruvka.c
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

    graph_adj_list output_graph(input_graph.num_vertices);

    bool completed = false;
    while (!completed) {

    }
    return output_graph;
}

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

std::vector<int> parent_arr_receive(int world, int rank, int parent_vector_size) {
    std::vector<int> r_vector(parent_vector_size);
    // Rank 0 --> receiving from all other processes
    if (rank == 0) {
        std::vector<std::vector<int>> tmp_vec(world - 1);
        for (int i = 1; i < world; i++) {
            // std::cout << "Rank 0 sending to " << i << std::endl;
            tmp_vec[i - 1].resize(parent_vector_size);
            MPI_Recv(&tmp_vec[i - 1][0], parent_vector_size, MPI_INT, i, MPI_PARENT_VEC_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // TODO Logic to determine which edge is the "Correct" edge
        for (auto x : tmp_vec) {
            for (auto y : x) {
                std::cout << y << ", ";
            }
            std::cout << std::endl;
        }
    } else { // Else we're receiving from 0 on initial broadcast
        MPI_Bcast(&r_vector[0], parent_vector_size, MPI_INT, 0, MPI_COMM_WORLD);
    }
    return r_vector;
}

void parent_arr_send(int world, int rank, std::vector<int> parent_vector) {
    if (rank == 0) {
        MPI_Bcast(&parent_vector[0], parent_vector.size(), MPI_INT, rank, MPI_COMM_WORLD);
    } else {
        std::cout << "Sending from rank " << rank << std::endl;
        MPI_Send(&parent_vector[0], parent_vector.size(), MPI_INT, 0, MPI_PARENT_VEC_TAG, MPI_COMM_WORLD);
    }
}


#endif

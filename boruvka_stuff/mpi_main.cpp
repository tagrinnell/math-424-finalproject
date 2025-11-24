
#include "include/graph_adj_list.hpp"
#include "include/boruvka.hpp"

#include <filesystem>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>

#include <mpi.h>

#define NUM_TRIALS 3

void mpi_test(int argc, char** argv);
graph_adj_list read_in_graph(std::string filepath);

int main(int argc, char** argv) {
    mpi_test(argc, argv);

    return 0;
}

void mpi_test(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    std::string graph_in = "graph_inputs/" + std::string(argv[argc - 1]);
    auto g = read_in_graph(graph_in);
    auto mpi = mpi_wrapper(g, argc, argv);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}


graph_adj_list read_in_graph(std::string filepath) {
    // std::cout << std::filesystem::current_path() << std::endl;
    graph_adj_list return_graph;

    std::ifstream f(filepath.c_str(), std::ios::in);
    int num_vertices = 0, num_edges = 0;

    if (!(f >> num_vertices >> num_edges)) {
        std::clog << "Error!! Graph failed" << std::endl;
        return return_graph;
    }

    return_graph.num_vertices = num_vertices;
    return_graph.init();

    int u, v, w;
    for (int i = 0; i < num_edges; i++) {
        f >> u >> v >> w;
        return_graph.add_edge(u, v, w);
    }

    return return_graph;
}
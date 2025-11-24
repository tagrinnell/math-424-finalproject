
#include "include/graph_adj_list.hpp"
#include "include/boruvka.hpp"

#include <filesystem>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#ifdef _OPENMP
    #include <omp.h>
#endif

void mpi_test(int argc, char** argv);
// bool assert_correctness(graph_adj_list a, graph_adj_list b);
graph_adj_list read_in_graph(std::string filepath);

int main(int argc, char** argv) {
    std::vector<std::string> test_inputs = {
        "graph10_25_50.txt",
        "graph10_40_20.txt",
        "graph1000_25000_10000.txt"
    };

    // std::cout << "Last arg " << argv[argc - 1] << std::endl;
    // std::string file = "graph_inputs/" + std::string(argv[argc - 1]);
    // auto graph_in = read_in_graph(file);
    // graph_in.to_string();
    mpi_test(argc, argv);

    return 0;
}

void mpi_test(int argc, char** argv) {
    // std::cout << "\n\nMPI TEST\n" << std::endl;
    graph_adj_list g(9);

    g.add_edge(0, 1, 17);
    g.add_edge(0, 2, 16);
    g.add_edge(0, 4, 2);
    g.add_edge(0, 5, 15);
    g.add_edge(4, 5, 3);
    g.add_edge(4, 1, 1);
    g.add_edge(5, 2, 2);
    g.add_edge(1, 2, 14);
    g.add_edge(1, 3, 14);
    g.add_edge(2, 3, 15);
    g.add_edge(2, 6, 4);
    g.add_edge(2, 7, 5);
    g.add_edge(3, 8, 7);
    g.add_edge(3, 7, 3);

    std::chrono::_V2::system_clock::time_point t1, t2;
    auto err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        std::cout << "MPI failed somehow" << std::endl;
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    t1 = std::chrono::high_resolution_clock::now();
    auto baseline = boruvka_mst(g);
    t2 = std::chrono::system_clock::now();

    auto mpi = mpi_wrapper(g, argc, argv);
    std::chrono::duration<double, std::milli> ms = t2 - t1;
    if (rank == 0) {
        std::cout << "Serial Implementation ran for " << ms.count() << std::endl;

        std::cout << "Asserting Correctness..." << std::flush;
        if (assert_correctness(baseline, mpi)) {
            std::cout << "Correct!" << std::endl;
        } else {
            std::cout << "Wrong :(" << std::endl;
        }

        std::cout << "Serial output: " << std::endl;
        baseline.to_string();
        std::cout << "MPI output: " << std::endl;
        mpi.to_string();
    }

    MPI_Finalize();
}

// bool assert_correctness(graph_adj_list a, graph_adj_list b) {
//     if (a.edge_list.size() != b.edge_list.size() && a.edge_list[0].size() != b.edge_list[0].size()) {
//         return false;
//     }
//     for (int i = 0; i < a.edge_list.size(); i++) {
//         for (int j = 0; j < a.edge_list[0].size(); j++) {
//             if (a.edge_list[i][j] != b.edge_list[i][j]) {
//                 return false;
//             }
//         }
//     }
//     return true;
// }

graph_adj_list read_in_graph(std::string filepath) {
    std::cout << std::filesystem::current_path() << std::endl;
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
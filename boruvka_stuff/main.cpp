
#include "graph.hpp"
#include "boruvka.hpp"

#include <chrono>
#include <iostream>
#ifdef _OPENMP
    #include <omp.h>
#endif

#ifdef _OPENMP
    #include <omp.h>
    graph boruvka_mst_openmp(graph input_graph, int num_threads);
#endif

void gfg_test();
void new_test();
graph boruvka_mst(graph input_graph);

int main() {
    // gfg_test();
    new_test();

    return 0;
}

void gfg_test() {
    graph g(4);

    g.add_edge(0, 1, 10);
    g.add_edge(0, 2, 6);
    g.add_edge(0, 3, 5);
    g.add_edge(1, 3, 15);
    g.add_edge(2, 3, 4);
    // g.add_edge(0, 1, 2);

    auto output = boruvka_mst(g);
    g.to_string();
    output.to_string();

}

void new_test() {
    std::cout << "\n\nNEW TEST\n" << std::endl;
    graph g(9);

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

    // g.add_edge(0, 1, 2);

    const auto t1 = std::chrono::high_resolution_clock::now();

    auto output = boruvka_mst(&g);
    const auto t2 = std::chrono::system_clock::now();
    std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << "Serial Implementation ran for " << ms.count() << std::endl;

    output.to_string();
    #ifdef _OPENMP
        std::cout << "\n\nOpenMP Trial" << std::endl;
        auto t3 = std::chrono::system_clock::now();
        output = boruvka_mst_openmp(&g, 16);
        auto t4 = std::chrono::system_clock::now();
        ms = t4 - t3;
        output.to_string();
        std::cout << "OpenMP Implementation ran for " << ms.count() << std::endl;
    #endif

}

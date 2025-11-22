
#include "boruvka.hpp"

// #define _MPI true

#include <vector>
#include <mpi.h>
#include <cmath>

// Parent array send / receive
void parent_arr_mpi(int size, int rank);
void edge_lists_mpi(int size, int rank);
void component_mpi(int size, int rank);
void print_vec(int rank, std::vector<int> vec);

// MPI functions
// std::vector<int> parent_arr_receive(int size, int rank, int parent_vector_size);
// void parent_arr_send(int size, int rank, std::vector<int> parent_vector);
// std::vector<int> edge_receive(int size, int rank, int dispatched_components);
// void edge_send(int size, int rank, std::vector<int> edges_to_send);
// std::unordered_map<int, bool> component_arr_receive(int size, int rank);
// void component_arr_send(int size, int rank, std::vector<int> component_list);

int main(int argc, char** argv) {
    auto err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        std::cout << "MPI ERROR" << std::endl;
        return -1;
    }

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    if (rank == 2) {
        MPI_Finalize();
        std::cout << "Rank 2 returning" << std::endl;
        return 0;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // srand(rank * size);
    // // std::cout << "Rank " << rank << ", world size " << size << std::endl;
    // if (rank == 0) {
    //     std::cout << "\n\t Parent Test" << std::endl;
    // }

    // parent_arr_mpi(size, rank);
    // MPI_Barrier(MPI_COMM_WORLD);
    // if (rank == 0) {
    //     std::cout << "\n\t Edge Test" << std::endl;
    // }
    // edge_lists_mpi(size, rank);
    // MPI_Barrier(MPI_COMM_WORLD);
    // if (rank == 0) {
    //     std::cout << "\n\t Component Test" << std::endl;
    // }
    // component_mpi(size, rank);
    // emulate_mpi(size, rank);

    MPI_Finalize();

    return 0;
}


// Parent array send / receive
void parent_arr_mpi(int size, int rank) {
    int num_indices = 20;
    if (rank == 0) {
        // std::cout << "\tPARENT TEST" << std::endl;
        // Generate random parent vector
        std::vector<int> parent_vec_send(num_indices);

        for (int i = 0; i < num_indices; i++) {
            parent_vec_send[i] = rand() * 1.0 / RAND_MAX * 20;
        }
        std::cout << "Root Rank sending: ";
        for (auto i : parent_vec_send) {
            std::cout << i << " ";
        }
        std::cout << std::endl;

        parent_arr_send(size, rank, parent_vec_send);


    } else {
        // Receive aa single array form ROOT_RANK
        auto rec_arr = parent_arr_receive(size, rank, num_indices);
        print_vec(rank, rec_arr);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

// TODO Edge lists send / receive
void edge_lists_mpi(int size, int rank) {
    if (rank == 0) {
        // std::cout << "\n\n\tEDGE TEST" << std::endl;
        // Generate random parent vector
        auto receive_edges = edge_receive(size, rank, size - 1);

        std::cout << "vec size " << receive_edges.size() << "  ";
        print_vec(rank, receive_edges);

    } else {
        // Generate a random vector of size 3x and send to root rank
        std::vector<int> random_vec(6 * rank);
        for (int i = 0; i < 6 * rank; i++) {
            random_vec[i] = 2 * i;
        }
        std::cout << "vec size " << 6 * rank << "  ";
        print_vec(rank, random_vec);
        edge_send(size, rank, random_vec);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

// TODO Component Array send / receive
void component_mpi(int size, int rank) {

    if (rank == 0) {
        std::vector<int> comp_list;
        for (int i = 0; i < 14; i++) {
            comp_list.push_back(i);
        }
        print_vec(rank, comp_list);
        std::cout << "  Expecting 14 / 4 = 3 with process 3  having 5 components" << std::endl;

        auto x = component_arr_scatter(size, rank, comp_list);
        std::cout << "Comp map from Rank " << rank << " ";
        for (auto y : x) {
            std::cout << y.first << " ";
        }
        std::cout << std::endl;
    } else {
        auto comp_map = component_arr_scatter(size, rank);
        std::cout << "Comp map from Rank " << rank << " ";
        for (auto x : comp_map) {
            std::cout << x.first << " ";
        }
        std::cout << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "\tComponent list size < world size" << std::endl;
        std::vector<int> comp_list;
        for (int i = 0; i < 3; i++) {
            comp_list.push_back(i);
        }
        print_vec(rank, comp_list);

        auto x = component_arr_scatter(size, rank, comp_list);
        std::cout << "Comp map from Rank " << rank << " ";
        for (auto y : x) {
            std::cout << y.first << " ";
        }
        std::cout << std::endl;
    } else {
        auto new_map = component_arr_scatter(size, rank);
        std::cout << "Comp map from Rank " << rank << " ";
        for (auto x : new_map) {
            std::cout << x.first << " ";
        }
        std::cout << std::endl;
    }
}

// TODO Do series of sends / receives like in the boruvka_mst_mpi function
void emulate_mpi(int size, int rank) {

}

void print_vec(int rank, std::vector<int> vec) {

    std::cout << "Rank " << rank << " printing vec: ";
    for (auto i : vec) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

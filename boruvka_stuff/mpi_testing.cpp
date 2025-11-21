
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
        return -1;
    }

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::cout << "Rank " << rank << ", world size " << size << std::endl;

    parent_arr_mpi(size, rank);

    MPI_Finalize();

    return 0;
}


// Parent array send / receive
void parent_arr_mpi(int size, int rank) {
    int num_indices = 20;
    if (rank == 0) {
        std::cout << "PARENT TEST" << std::endl;
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

}

// TODO Component Array send / receive
void component_mpi(int size, int rank) {

}

// TODO Do series of sends / receives like in the boruvka_mst_mpi function
void emulate_mpi() {

}

void print_vec(int rank, std::vector<int> vec) {

    std::cout << "Rank " << rank << " printing vec: ";
    for (auto i : vec) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

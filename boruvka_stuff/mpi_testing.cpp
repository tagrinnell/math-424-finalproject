
#include "boruvka.hpp"
#include <mpi.h>

int main(int argc, char** argv) {
    auto err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        return -1;
    }
    std::vector<int> x = {0, 1, 2, 3, 4, 5};
    std::vector<int> y = {13, 14, 15, 12, 11, 0};
    std::vector<int> z = {-1, -2, -3, 5, 2, 1};

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::cout << "Rank: " << rank << std::endl;
    std::cout << "Size: " << size << std::endl;

    // Test broadcasting from 0 -> 1 / 2
    if (rank == 0) {
        std::cout << "Sending from rank 0" << std::endl;
        // MPI_Send(&x[0], x.size(), MPI_INT, 1, 0, MPI_COMM_WORLD);

        parent_arr_send(size, rank, x);
    } else {
        // std::vector<int> buff(x.size());
        // MPI_Recv(&buff[0], x.size(), MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // for (auto i : buff) {
        //     std::cout << i << " ";
        // }
        // std::cout << "\nReceive First " << std::endl;
        // MPI_Send(&x[0], sizeof(x), MPI_INT, 1, MPI_ANY_TAG, MPI_COMM_WORLD);
        auto d = parent_arr_receive(size, rank, x.size());
        for (auto i : d) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Test Gathering from
    if (rank == 0) {
        std::cout << "Testing 'Gathering'" << std::endl;
        parent_arr_receive(size, rank, x.size());
    } else {
        if (rank == 1) {
            parent_arr_send(size, rank, y);
        } else if (rank == 2) {
            parent_arr_send(size, rank, z);
        }
    }

    MPI_Finalize();

    return 0;
}

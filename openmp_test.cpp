
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

int main() {

#ifdef _OPENMP
    #pragma omp parallel
    {
        std::cout << "Outputting from thread " << omp_get_thread_num() << std::endl;
    }
#endif

#ifndef _OPENMP
    std::cout << "Failed to set _OPENMP Flag" << std::endl;
#endif

    return 0;
}
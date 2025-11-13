#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <chrono>

#define NUM_ITERATIONS 100

long double vector_benchmark() {

    long double total_sum = 0.0;
    for (int i = 0; i < NUM_ITERATIONS; i++) {

        auto t1 = std::chrono::high_resolution_clock::now();

        std::vector<int> x = {0, 1, 2};
        int y = x[1];

        auto t2 = std::chrono::high_resolution_clock::now();
        auto d = std::chrono::duration<double, std::milli>(t2 - t1);
        std::cout << d << std::endl;
    }

    return total_sum / (1.0 * NUM_ITERATIONS);
}

long double tuple_benchmark() {

    long double total_sum = 0.0;
    for (int i = 0; i < NUM_ITERATIONS; i++) {

        auto t1 = std::chrono::high_resolution_clock::now();

        std::tuple<int, int, int> x = std::make_tuple(0, 1, 2);
        int y = std::get<1>(x);

        auto t2 = std::chrono::high_resolution_clock::now();
        total_sum += (t2 - t1).count();
    }

    return total_sum / (1.0 * NUM_ITERATIONS);
}

int main() {

    auto tuple_res = tuple_benchmark();
    auto vec_res = vector_benchmark();

    std::cout << "Tuple ran " << NUM_ITERATIONS << " Iterations in " << tuple_res << "ms" << std::endl;
    std::cout << "Vector ran " << NUM_ITERATIONS << " Iterations in " << vec_res << "ms" << std::endl;

    return 0;
}

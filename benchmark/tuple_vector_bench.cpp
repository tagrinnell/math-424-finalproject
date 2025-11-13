#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <chrono>

#define NUM_ITERATIONS 1000000

double vector_benchmark() {

    float total_sum = 0.0;
    for (int i = 0; i < NUM_ITERATIONS; i++) {

        auto t1 = std::chrono::system_clock::now();

        std::vector<int> x = {0, 1, 2};
        int y = x[1];

        auto t2 = std::chrono::system_clock::now();
        total_sum += (t2 - t1).count();
    }

    return total_sum / (1.0 * NUM_ITERATIONS);
}

double tuple_benchmark() {

    float total_sum = 0.0;
    for (int i = 0; i < NUM_ITERATIONS; i++) {

        auto t1 = std::chrono::system_clock::now();

        std::tuple<int, int, int> x = std::make_tuple(0, 1, 2);
        int y = std::get<1>(x);

        auto t2 = std::chrono::system_clock::now();
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

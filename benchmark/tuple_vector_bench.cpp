#include <iostream>
#include <vector>
#include <thread>
#include <tuple>
#include <cmath>
#include <chrono>

#define NUM_ITERATIONS 1000

long double vector_benchmark() {
    std::chrono::_V2::system_clock::time_point t1, t2;
    double total_sum = 0.0;
    for (int i = 0; i < NUM_ITERATIONS; i++) {

        t1 = std::chrono::high_resolution_clock::now();

        std::vector<int> x = {0, 1, 2};
        int y = x[1];

        t2 = std::chrono::high_resolution_clock::now();
        // double x =
        total_sum += ((const std::chrono::duration<double, std::micro>)(t2 - t1)).count();
        std::cout << "Iteration " << i << " ran in " << ((const std::chrono::duration<double, std::milli>)(t2 - t1)).count() << " micro-seconds" << std::endl;
    }

    // const std::chrono::duration<double, std::milli> duration = t2 - t1;
    // std::cout << "Outputting " << duration.count() << "milliseconds" << std::endl;
    return total_sum / (1.0 * NUM_ITERATIONS);
}

long double tuple_benchmark() {

    long double total_sum = 0.0;
    for (int i = 0; i < NUM_ITERATIONS; i++) {

        auto t1 = std::chrono::high_resolution_clock::now();

        std::tuple<int, int, int> x = std::make_tuple(0, 1, 2);
        int y = std::get<1>(x);

        auto t2 = std::chrono::high_resolution_clock::now();
        total_sum += ((const std::chrono::duration<double, std::milli>)(t2 - t1)).count();
    }

    return total_sum / (1.0 * NUM_ITERATIONS);
}

int main() {
    auto t1 = std::chrono::high_resolution_clock::now();
    std::this_thread::sleep_for(std::chrono::seconds(1));
    auto t2 = std::chrono::high_resolution_clock::now();

    const auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "t2 - t1 = " << (int_ms).count() << " milliseconds" << std::endl;

    // auto tuple_res = tuple_benchmark();
    // auto vec_res = vector_benchmark();

    // std::cout << "Tuple ran " << NUM_ITERATIONS << " Iterations averaging " << tuple_res << "ms" << std::endl;
    // std::cout << "Vector ran " << NUM_ITERATIONS << " Iterations averaging " << vec_res << "ms" << std::endl;

    return 0;
}

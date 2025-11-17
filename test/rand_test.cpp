#include <vector>
#include <iostream>

int main() {
    std::vector<int> x(6);
    x[2] = 3;
    x[5] = 2;

    std::cout << x[2] << " " << x[5] << std::endl;
    for (int i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "]: " << x[i] << std::endl;
    }
}
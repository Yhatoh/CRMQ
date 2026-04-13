#include <vector>
#include <iostream>
#include <random>

#include "huff_vector.hpp"

int main(void) {

    //std::vector<int64_t> data = {2, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 10};

    std::vector<int> vec(1e6);
    std::random_device rd;  // Seed
    std::mt19937 gen(rd()); // Mersenne Twister RNG
    std::uniform_int_distribution<> dist(0, 100000);

    for (auto& val : vec) {
        val = dist(gen);
    }

    huff_vector<int> hv(vec);

    for(uint64_t i = 0; i < vec.size(); ++i)
        assert(vec[i] == hv[i]);    
    //std::cout << "expected: " << val[i] << " obtained: " << hv[i] << std::endl;

    return 0;
}

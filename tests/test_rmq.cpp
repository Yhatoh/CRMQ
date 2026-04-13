#include <iostream>
#include <string>
#include <vector>

#include "crmq.hpp"

int main() {

    std::vector<int> data = {3, 2, 4, 1, 7, 6, 8, 5, 9, 0, 12, 11, 13, 10, 14};

    CRMQ<int> crmq(data);

    std::cout << "rmq(5, 10): " << crmq.query(5, 10) << std::endl;

    return 0;
}
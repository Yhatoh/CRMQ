// std includes
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>
#include <fstream>
#include <filesystem>
#include <numeric>
#include <fstream>
#include <sstream>
#include <stack>
#include <map>

// local includes
#include "../libsais/include/libsais64.h"
#include "../libsais/include/libsais.h"
#include "crmq.hpp"

template<size_t threshold = 8>
void measure_bpe(const std::string &dataset, std::vector<int64_t> &data,
                                std::ofstream &csv_output) {
    const size_t n = data.size();

    auto start = std::chrono::high_resolution_clock::now();

    CRMQ<int64_t, threshold> crmq(data);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;

    const double time = duration.count();
    
    std::map<std::string, uint64_t> sizes = crmq.size_per_component();

    for(auto &comp : sizes) {
        csv_output << dataset << "," << time << "," << threshold << ","
                << comp.first << "," << double(comp.second) / double(n) << std::endl;
    }

    csv_output << dataset << "," << time << "," << threshold << "," 
                << "whole" << "," << double(crmq.size()) / double(n) << std::endl;

    std::cout << "Finished measuring the bpe for " << dataset << std::endl;
}

int main(int argc, char *argv[]) {
    if(argc < 2) {
        std::cerr << "Please provide the dataset directory and the results path" << std::endl;
        std::cerr << "Usage: ./measure_bpe_lcp datasets results_path" << std::endl;
        return -1;
    }

    const std::string ds_path = argv[1];
    const std::string out_dir = argv[2];

    std::ofstream csv_output(out_dir);

    if (!csv_output.is_open()) 
        throw std::runtime_error("Unable to open the file");

	csv_output << "dataset," << "time," << "min_nodes," << "component," << "bpe" << std::endl;
 
    for (const auto& entry : std::filesystem::directory_iterator(ds_path)) {
        const std::string path = entry.path();
        const std::string dataset = entry.path().filename();

        std::cout << "reading " << dataset << std::endl;

        if (std::filesystem::is_regular_file(path)) {

            std::ifstream file(path, std::ios::in | std::ios::binary);
    
            if(!file) throw std::runtime_error("Unable to open the file");
        
            std::string s((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

            const size_t n = s.size();
            int64_t *lcp_ptr;

            // build the suffix array and the permuted lcp array
            uint8_t *p = reinterpret_cast<uint8_t*>(s.data());
            int64_t *sa = new int64_t[n];
            int64_t *plcp = new int64_t[n];
            lcp_ptr = new int64_t[n];

            if(libsais64(p, sa, n, 0, NULL) != 0) throw std::runtime_error("SA construction failed");

            if(libsais64_plcp(p, sa, plcp, n) != 0) throw std::runtime_error("PLCP array construction failed");
            
            if(libsais64_lcp(plcp, sa, lcp_ptr, n) != 0) throw std::runtime_error("LCP array construction failed");

            std::cout << "Finished building the lcp array for " << dataset << std::endl;

            std::vector<int64_t> lcp;
            lcp.assign(lcp_ptr, lcp_ptr + n);

            measure_bpe<32>(dataset, lcp, csv_output);
            measure_bpe<48>(dataset, lcp, csv_output);
            measure_bpe<64>(dataset, lcp, csv_output);
		}
    }

    csv_output.close();

    return 0;
}

#include <gtest/gtest.h>
#include <random>

#include <sdsl/bit_vectors.hpp>

#include "scst.hpp"

class CCTTest : public ::testing::Test {
protected:

    static sdsl::bit_vector bp;

    static scst<> cct;

    static constexpr int seed = 42;

    static constexpr uint64_t size = 1e7;

    static sdsl::bit_vector build_bp_sequence(std::vector<int> &data) {
        data.push_back(std::numeric_limits<int>::min());
        
        const size_t n = data.size();
        sdsl::bit_vector bp(2 * n, 0);

        std::stack<size_t> st;
        size_t pos = 2 * n;
        for(int64_t i = n - 1; i >= 0; i--) {
            while (!st.empty() && data[st.top()] >= data[i]) {
                st.pop();
                bp[--pos] = 1;
            }
            st.push(i);
            bp[--pos] = 0;
        }
        while (!st.empty()) {
            st.pop();
            bp[--pos] = 1;
        }

        data.pop_back();

        return bp;
    }

    static void init_bp() {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<int> dis(1, size);

        std::vector<int> data;

        for (size_t i = 0; i < size; ++i)
            data.push_back(dis(gen));

        bp = build_bp_sequence(data);
        cct = scst<>(bp, 2);
    }

    static void SetUpTestSuite() {
        if(!bp.size()) {
            init_bp();  
        }
    }
};

sdsl::bit_vector CCTTest::bp;
scst<> CCTTest::cct;

TEST_F(CCTTest, InnodeTest) {
    sdsl::select_support_mcl<0> select0(&CCTTest::bp);
    for(uint64_t i = 1; i <= CCTTest::size; ++i) {
        uint64_t expected_innode = select0(i + 1);
		uint64_t new_computed_innode = CCTTest::cct.select0(i + 1);
        ASSERT_EQ(new_computed_innode, expected_innode) << " Computed: " 
                << new_computed_innode << ", Expected: " << expected_innode;
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

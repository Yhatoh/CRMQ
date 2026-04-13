#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include <utility>
#include <random>
#include <tuple>

#include "crmq.hpp"

using query_type = std::pair<uint64_t, uint64_t>;

template<typename K>
inline std::pair<K, uint64_t> real_minimum(const std::vector<K> &data, const uint64_t lo, const uint64_t hi) {
    K min = data[lo];
    uint64_t idx = lo;

    for(uint64_t i = lo + 1; i <= hi; ++i)
        if(data[i] < min) {
            min = data[i];
            idx = i;
        }

    return std::make_pair(min, idx);
}

typedef ::testing::Types<CRMQ<int32_t, 8>, CRMQ<int32_t, 16>, CRMQ<int32_t, 32>> CRMQTypes;

template <typename T>
class RMQTest : public ::testing::Test {
protected:

    static std::vector<int32_t> data;

    static std::vector<query_type> queries;

    static constexpr int seed = 42;

    static void init_data(uint64_t size) {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<int32_t> dis(1, size);

        for(uint64_t i = 0; i < size; ++i)
            data.push_back(dis(gen));
    }

    static void init_queries(uint64_t num_queries, uint64_t size) {
        const std::vector<uint64_t> lenghts = {200, 600, 1200, 
                                                6000, 8000, 10000};
        queries.reserve(num_queries * lenghts.size());
        std::mt19937 gen(seed);

        for(const auto &length : lenghts) {
            std::uniform_int_distribution<uint64_t> length_dis(1, size - length + 1);

            for(auto q = 0; q < num_queries; ++q) {
                const uint64_t start = length_dis(gen);
                const uint64_t end = start + length;
                queries.emplace_back(start, end);
            }
        }
    }

    static void SetUpTestSuite() {
        if(queries.empty()) {
            uint64_t size = 1e7;
            init_data(size);
            init_queries(10000, size);
        }
    }
};

template<typename T>
std::vector<int32_t> RMQTest<T>::data;

template<typename T>
std::vector<query_type> RMQTest<T>::queries;

template <typename T>
class CRMQTest : public RMQTest<T> {
protected:
    static void SetUpTestSuite() {
        RMQTest<T>::SetUpTestSuite();
    }
};

TYPED_TEST_SUITE(CRMQTest, CRMQTypes);

TYPED_TEST(CRMQTest, TestRMQs) {
    TypeParam rmq_ds(CRMQTest<TypeParam>::data);
    for(const auto &q : CRMQTest<TypeParam>::queries) {
        const auto [expected_min, expected_min_pos] = real_minimum(CRMQTest<TypeParam>::data, q.first, q.second);
        const auto computed_min_pos = rmq_ds.query(q.first, q.second);
        ASSERT_EQ(expected_min_pos, computed_min_pos) << " Query i = " << q.first << ", j = " << q.second;
        ASSERT_EQ(expected_min, CRMQTest<TypeParam>::data[computed_min_pos]) << " Query i = " << q.first << ", j = " << q.second;
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

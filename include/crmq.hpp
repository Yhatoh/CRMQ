#ifndef __CRMQ__
#define __CRMQ__
//  c++ includes
#include <cstdint>
#include <climits>
#include <string>
#include <vector>
#include <map>

// sdsl includes
#include <sdsl/int_vector.hpp>

// local includes
#include "scst.hpp"
#include "sdsl/io.hpp"
#include "sparse_table.hpp"

template<typename T, size_t threshold = 8>
class CRMQ {
	scst<> cct;
	size_t n;

	sdsl::int_vector<> min_data_samples;
	sdsl::int_vector<> pos_min_data_samples;

	sparse_table<sdsl::int_vector<>> top_level_st;

	// only for paper purposes
	uint64_t ans_sparse_table;

	public:

	uint64_t get_sparse_table() { return ans_sparse_table; }
	void reset_sparse_table() { ans_sparse_table = 0; }

	uint64_t get_recursive_calls() { return cct.get_recursive_calls(); }
	void reset_recursive_calls() { cct.reset_recursive_calls(); }

	CRMQ() = default;

	explicit CRMQ(sdsl::bit_vector bp) : n(bp.size() / 2) {
	}

	explicit CRMQ(std::vector<T> &data) : n(data.size()) {
		if(n == 0) [[unlikely]]
			return;

		sdsl::bit_vector bp = build_bp_sequence(data);

		auto aux = scst<>(bp, threshold);

		cct = aux;

		sample_data(data);
		top_level_st = sparse_table<sdsl::int_vector<>>(&min_data_samples);
	}

	[[nodiscard]] inline uint64_t query(const uint64_t i, const uint64_t j) {
		const uint64_t DataSamples = 2048;
		if constexpr (DataSamples > 0) {
			const size_t ds_i = i / DataSamples;
			const size_t ds_j = j / DataSamples;

			const size_t min_sample = top_level_st.query(ds_i, ds_j);
			const size_t min_idx = pos_min_data_samples[min_sample] + min_sample * DataSamples;

			if(i <= min_idx && min_idx <= j) {
			  ans_sparse_table++;
                          return min_idx;
                        }
		}

		uint64_t ii = cct.select0(i + 1);
		uint64_t jj = cct.select0(j + 1);

		assert(ii >= 0 && ii < 2 * n + 2);
		assert(jj >= 0 && jj < 2 * n + 2);

		uint64_t c = 0;

		const uint64_t lca = cct.lca(ii, jj, 0, c);

		return cct.inorder(lca, c);
	}

	inline uint64_t size() {
		return cct.size_in_bits() + (sizeof(n) * CHAR_BIT)
			+ top_level_st.size() +
			(sdsl::size_in_bytes(min_data_samples) +
			 sdsl::size_in_bytes(pos_min_data_samples)) * CHAR_BIT;
	}

	inline std::map<std::string, uint64_t> size_per_component() {
		auto ret = cct.size_per_component();
		ret["12. sparsetable sparse table"] = top_level_st.size();
		ret["13. sparsetable min data samples"] = sdsl::size_in_bytes(min_data_samples) * CHAR_BIT;
		ret["14. sparsetable pos min data samples"] = sdsl::size_in_bytes(pos_min_data_samples) * CHAR_BIT;
		return ret;
	}

	private: 

	sdsl::bit_vector build_bp_sequence(std::vector<T> &data) const {
		// this encloses everything into an "()"
		// todo: can we avoid to change the array? this does not
		// work if the input already contains such a value
		//data.push_back(std::numeric_limits<T>::min());

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

		//data.pop_back();

		return bp;
	}

	void sample_data(const std::vector<T> &data) {
		const uint64_t DataSamples = 2048;
		auto samples = n / DataSamples + ((n % DataSamples) != 0);

		std::vector<T> tmp_min_data_samples(samples);
		min_data_samples = sdsl::int_vector<>(samples);
		pos_min_data_samples = sdsl::int_vector<>(samples);

		for(auto i = 0; i < samples; ++i) {
			auto min_idx = i * DataSamples;
			for(auto j = min_idx; j < (i+1) * DataSamples && j < n; ++j) {
				if(data[j] < data[min_idx]) min_idx = j; // notice: leftmost minimum here
			}
			tmp_min_data_samples[i] = data[min_idx];
			pos_min_data_samples[i] = min_idx - i * DataSamples;
		}

		// since data can contain arbitrarily large values
		// we remap its content into [1, n / DataSamples].
		// Notice that this does not change the result of any query

		std::vector<T> sorted_samples = tmp_min_data_samples;
		std::sort(sorted_samples.begin(), sorted_samples.end());

		std::unordered_map<T, size_t> rank_map;
		size_t rank = 0;

		for(const auto &sample : sorted_samples) {
			if(rank_map.find(sample) == rank_map.end())
				rank_map[sample] = rank++;
		}

		for(auto i = 0; i < tmp_min_data_samples.size(); ++i)
			min_data_samples[i] = rank_map[tmp_min_data_samples[i]];

		sdsl::util::bit_compress(min_data_samples);
		sdsl::util::bit_compress(pos_min_data_samples);
	}
};

#endif

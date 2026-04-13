#ifndef __CBP__
#define __CBP__
//  c++ includes
#include <cstdint>
#include <vector>

// sdsl includes
#include <sdsl/bp_support_sada.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>

// local includes
#include "../libsais/include/libsais64.h"
#include "../libsais/include/libsais.h"
#include "union_find.hpp"

class cbp {
	sdsl::bit_vector cbv;
	// onlye useful for debugging purposes
	sdsl::bp_support_sada<> support_cbv;

	sdsl::int_vector<> p;

	// bit_vector s
	sdsl::sd_vector<> s;
	sdsl::rank_support_sd<> rank_s;
	sdsl::select_support_sd<> select_s;

	// bit vector r
	sdsl::sd_vector<> r;
	sdsl::rank_support_sd<> rank_r;
	sdsl::select_support_sd<> select_r;

	public:

	uint64_t size() { return cbv.size(); }
	uint64_t max_trees() { return p.size(); }
	uint64_t nodes() { return cbv.size() / 2; }

	cbp(sdsl::bit_vector &bv, const std::uint32_t threshold = 8) {
		sdsl::bp_support_sada<> tree_support = sdsl::bp_support_sada<>(&bv);

		uint8_t* text = new uint8_t[bv.size()];
		for(uint64_t i = 0; i < bv.size(); i++) {
			text[i] = bv[i];
		}

		int64_t *csa = new int64_t[bv.size()];
		int64_t *plcp = new int64_t[bv.size()];
		int64_t *lcp = new int64_t[bv.size()];

		if(libsais64(text, csa, bv.size(), 0, NULL) != 0)
			throw std::runtime_error("SA construction failed");
		if(libsais64_plcp(text, csa, plcp, bv.size()) != 0)
			throw std::runtime_error("PLCP array construction failed");
		if(libsais64_lcp(plcp, csa, lcp, bv.size()) != 0)
			throw std::runtime_error("LCP array construction failed");

		union_find idems_tree(bv.size());
		for(uint64_t i = 1; i < bv.size(); i++) {
			// only considering suffix starting with (
			if(bv[csa[i]]) {

				uint64_t curr_open = csa[i];
				uint64_t curr_close = tree_support.find_close(curr_open);

				uint64_t prev_open = csa[i - 1];

				// ignoring leaves
				if(curr_close - curr_open + 1 <= threshold) continue;

				if(lcp[i] >= curr_close - curr_open + 1) {
					idems_tree.union_set(curr_open, prev_open);
				}
			}
		}

		// clean, is completly useless
		free(csa);
		free(plcp);
		free(lcp);

		uint64_t bit = 0;
		std::vector< uint64_t > compr_bv;
		std::vector< uint64_t > p_aux;

		std::vector< uint64_t > s_bv;
		std::vector< uint64_t > r_bv;

		uint64_t amount_of_bits_removed = 0;
		// this is to compute the correct position
		std::vector< uint64_t > prefix_help(bv.size(), 0);

		for(uint64_t i = 0; i < bv.size(); i++) {
			if(bv[i]) {
				compr_bv.push_back(bit++);

				uint64_t leader = idems_tree.find_set(i);
				// cutting the subtree
				if(leader != i) {
					r_bv.push_back(bit - 1);
					s_bv.push_back(i);
					p_aux.push_back(leader - prefix_help[leader - 1]);
					bit++;

					uint64_t next_i = tree_support.find_close(i);
					for(uint64_t pfh = i; pfh <= next_i; pfh++) {
						prefix_help[pfh] = prefix_help[pfh - 1] + 1;
					}
					prefix_help[next_i] -= 2;
					i = next_i;
					s_bv.push_back(i);
				} else {
					if(i > 0) prefix_help[i] = prefix_help[i - 1];
				}
			} else {
				prefix_help[i] = prefix_help[i - 1];
				bit++;
			}
		}

		prefix_help.clear();
		idems_tree.clear();

		// creation of the cbv
		cbv = sdsl::bit_vector(bit, 0);
		for(const auto& i : compr_bv) cbv[i] = 1;
		support_cbv = sdsl::bp_support_sada<>(&cbv);

		// creation of int_vector for positions
		p = sdsl::int_vector<>(p_aux.size(), 0);
		for(uint64_t i = 0; i < p_aux.size(); i++) {
			p[i] = p_aux[i];
		}

		sdsl::util::bit_compress(p);

		s = sdsl::sd_vector<>(s_bv.begin(), s_bv.end());
		sdsl::util::init_support(rank_s, &s);
		sdsl::util::init_support(select_s, &s);

		r = sdsl::sd_vector<>(r_bv.begin(), r_bv.end());
		sdsl::util::init_support(rank_r, &r);
		sdsl::util::init_support(select_r, &r);
	}

	void decompress(sdsl::bit_vector &bp) {
		std::vector< uint64_t > ones;
		uint64_t bit = 0;
		std::stack< std::pair< int, int > > rec;
		for(uint64_t i = 0; i < cbv.size(); i++) {
			if(cbv[i]) {
				if(r[i]) {
					uint64_t r_i = rank_r(i);
					rec.push({i, support_cbv.find_close(p[r_i])});
					i = p[r_i];
				}
				ones.push_back(bit++);
			} else {
				if(!rec.empty() && i >= rec.top().second) {
					i = rec.top().first + 1;
					rec.pop();
				}
				bit++;
			}
		}
		bp = sdsl::bit_vector(bit, 0);
		for(auto &i : ones) bp[i] = 1;
	}

};

#endif

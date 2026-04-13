#ifndef __SCST__
#define __SCST__
//  c++ includes
#include <cstdint>
#include <climits>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <set>

// sdsl includes
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/vlc_vector.hpp>

// pasta includes
#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/find_l2_flat_with.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>

// local includes
#include "my_bp_support.hpp"
#include "../libsais/include/libsais64.h"
#include "bp_utils.hpp"
#include "sux/bits/EliasFano.hpp"

template< uint64_t select_bs = 64, uint64_t exc_sample_bs = 2048 >
class scst {
	private:

		pasta::BitVector cbv;
		pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES> support_cbv;

		sdsl::int_vector<> p_s;

		uint64_t first_pointer_s;
		uint64_t first_pointer_r;

		sdsl::int_vector<> exc_samples;
		sdsl::int_vector<> min_exc_samples;
		sdsl::int_vector<> pos_min_exc_samples;

		sux::bits::EliasFano<> r;
		sux::bits::EliasFano<> s;

		uint64_t innode_pointers_ones;
		sux::bits::EliasFano<> innode_pointers;

		// only for paper purposes
		uint64_t recursive_calls;

		inline uint64_t phi_1(uint64_t u) {
			uint64_t r_value = s.rank(u);
			if(r_value & 1) {
				return r.select((r_value + 1) / 2 - 1);
			}

			if(r_value == 0) return u;
			else
				return r.select(r_value / 2 - 1) + u - s.select(r_value - 1);
		}

		inline uint64_t phi(uint64_t v) {
			uint64_t r_value = r.rank(v);
			if(r_value == 0) return v;
			else {
				if(r.rank(v + 1) > r_value)
					return s.select(2 * r_value - 1) + v - r.select(r_value - 1) - 1;
				else
					return s.select(2 * r_value - 1) + v - r.select(r_value - 1);
			}
		}

	public:

		void reset_recursive_calls() { recursive_calls = 0; }
		uint64_t get_recursive_calls() { return recursive_calls; }

		uint64_t max_trees() { return p_s.size(); }
		uint64_t nodes() { return cbv.size() / 2; }

		scst() = default;

		scst(sdsl::bit_vector &bv, const std::uint32_t threshold = 2) {
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

			free(plcp);
			free(text);

			std::vector< uint64_t > idems_tree(bv.size(), 0);
			for(uint64_t i = 0; i < bv.size(); i++)
				idems_tree[i] = i;
			uint64_t start_idem = -1;
			uint64_t min_idx = 0;
			for(uint64_t i = 1; i < bv.size(); i++) {
				// only considering suffix starting with (
				if(bv[csa[i]]) {

					uint64_t curr_open = csa[i];

					uint64_t curr_close = tree_support.find_close(curr_open);
					uint64_t enclose = (curr_close == bv.size() - 1) ? bv.size() - 1 : tree_support.fwd_excess(curr_close, -1);
					uint64_t subtree_size = enclose - curr_open;

					if((lcp[i] >= subtree_size) &&
							(subtree_size > threshold + 2 + 1 + 1)) {
						if(start_idem == -1) {
							start_idem = i - 1;
							min_idx = std::min(csa[start_idem], csa[i]);
						} else {
							min_idx = std::min(min_idx, (uint64_t) csa[i]);
						}
					} else if(start_idem != 0) { // i found a set of identical subtrees
						for(uint64_t l = start_idem; l < i; l++) {
							idems_tree[csa[l]] = min_idx;
						}
						start_idem = -1;
					}

				}
			}

			// clean, is completly useless
			free(csa);
			free(lcp);

			uint64_t bit = 0;
			std::vector< uint64_t > compr_bv;
			std::vector< uint64_t > p_aux;
			std::vector< uint64_t > p_s_aux;
			std::vector< uint64_t > s_bv;
			std::vector< uint64_t > r_bv;

			uint64_t amount_of_bits_removed = 0;
			// this is to compute the correct position
			std::vector< uint64_t > prefix_help(bv.size(), 0);

			uint64_t innode = 0;
			std::vector< uint64_t > pb_innodes;
			pb_innodes.push_back(0);

			for(uint64_t i = 0; i < bv.size(); i++) {
				if(bv[i]) {
					compr_bv.push_back(bit++);

					uint64_t leader = idems_tree[i];
					// cutting the subtree
					if(leader != i) {
						r_bv.push_back(bit);
						s_bv.push_back(i);
						p_aux.push_back(leader - prefix_help[leader - 1]);
						p_s_aux.push_back(leader);
						bit++;

						uint64_t curr_close = tree_support.find_close(i);
						uint64_t enclose =
							(curr_close == bv.size() - 1) ?
							bv.size() - 1 :
							tree_support.fwd_excess(curr_close, -1);
						uint64_t subtree_size = enclose - i;
						innode += subtree_size / 2;
						uint64_t next_i = i + subtree_size - 1;

						for(uint64_t pfh = i; pfh <= next_i; pfh++) {
							prefix_help[pfh] = prefix_help[pfh - 1] + 1;
						}
						prefix_help[next_i] -= 2;
						i = next_i;
						pb_innodes.push_back(innode);
						s_bv.push_back(i);
					} else {
						if(i > 0) prefix_help[i] = prefix_help[i - 1];
					}
				} else {
					innode++;
					if(i > 0) prefix_help[i] = prefix_help[i - 1];
					bit++;
				}
			}

			prefix_help.clear();
			idems_tree.clear();

			cbv = pasta::BitVector(bit, 0);
			for(const auto& i : compr_bv) cbv[i] = 1;
			support_cbv = pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES>(cbv);

			std::vector< int64_t > excess_array(bit, 0);
			int64_t exc = 0;
			for(uint64_t i = 0; i < bit; i++) {
				exc += (cbv[i] ? 1 : -1);
				excess_array[i] = exc;
			}
			sample_excess_array(excess_array);

			// creation of int_vector for positions
			p_s = sdsl::int_vector<>(p_s_aux.size(), 0);
			for(uint64_t i = 0; i < p_s_aux.size(); i++) {
				p_s[i] = p_s_aux[i];
			}
			sdsl::util::bit_compress(p_s);

			s = sux::bits::EliasFano<>(s_bv, bv.size());
			r = sux::bits::EliasFano<>(r_bv, cbv.size());
			innode_pointers_ones = pb_innodes.size();
			innode_pointers = sux::bits::EliasFano<>(pb_innodes, pb_innodes.back() + 1);

			// todo: don't use select for this
			first_pointer_s = s.select(0);
			first_pointer_r = r.select(0) - 1; // select_r(0) should never be < 1
		}

		uint64_t lca(uint64_t u, uint64_t v) {
			if(u > v) std::swap(u, v);

			uint64_t r_u = s.rank(u);
			uint64_t r_v = s.rank(v);

			// they are not pointers
			// or just one
			if((r_u & 1) + (r_v & 1) < 2) {
				uint64_t u_ = phi_1(u);
				uint64_t v_ = phi_1(v);
				if(v_ - u_ + 1 <= exc_sample_bs) {
					const int64_t exc_i = u_ - 2 * support_cbv.rank0(u_);
					auto min_exc = min_excess_scan(cbv, exc_i, u_, v_);
					const uint64_t lca_idx = phi(min_exc.first);
					return lca_idx;
				}
				auto [lca_idx_c, min_got] = min_excess(u_, v_);
				const uint64_t lca_idx = phi(lca_idx_c);

				return lca_idx;
			}

			// both pointers but not the same
			if(r_u != r_v) {
				uint64_t u_ = phi_1(u);
				uint64_t v_ = phi_1(v);
				if(v_ - u_ + 1 <= exc_sample_bs) {
					const int64_t exc_i = u_ - 2 * support_cbv.rank0(u_);
					auto min_exc = min_excess_scan(cbv, exc_i, u_, v_);
					const uint64_t lca_idx = phi(min_exc.first);
					return lca_idx;
				}
				auto [lca_idx_c, min_got] = min_excess(u_, v_);
				const uint64_t lca_idx = phi(lca_idx_c);

				return lca_idx;
			}

			uint64_t p = r_u;
			uint64_t m = s.select(p - 1);
			uint64_t rank_tree = ((p + 1) / 2) - 1;
			uint64_t u_ = p_s[rank_tree] + u - m;
			uint64_t v_ = p_s[rank_tree] + v - m;
			return lca(u_, v_) - p_s[rank_tree] + m;
		}

		uint64_t lca(uint64_t u, uint64_t v, const uint64_t b, uint64_t &c) {
			if(u > v) [[unlikely]] std::swap(u, v);

			//assert(u >= 0 && u < s.size());
			//assert(v >= 0 && v < s.size());

			const uint64_t r_u = s.rank(u);
			const uint64_t r_v = s.rank(v);

			// todo: we can postpone the rank operations to the inorder function
			// so as to avoid them we they are not necessary

			// they are not pointers
			// or just one
			if((r_u & 1) + (r_v & 1) < 2) {
				uint64_t u_ = 0;
				if(r_u & 1) {
					u_ = r.select((r_u + 1) / 2 - 1);
				} else if(r_u == 0) u_ = u;
				else u_ = r.select(r_u / 2 - 1) + u - s.select(r_u - 1);

				uint64_t v_ = 0;
				if(r_v & 1) {
					v_ = r.select((r_v + 1) / 2 - 1);
				} else if(r_v == 0) v_ = v;
				else v_ = r.select(r_v / 2 - 1) + v - s.select(r_v - 1);

				if(v_ - u_ + 1 <= exc_sample_bs) {
					const int64_t exc_i = u_ - 2 * support_cbv.rank0(u_);
					auto min_exc = min_excess_scan(cbv, exc_i, u_, v_);
					c = support_cbv.rank0(min_exc.first) - support_cbv.rank0(b);
					const uint64_t lca_idx = phi(min_exc.first);
					return lca_idx;
				}
				auto [lca_idx_c, min_got] = min_excess(u_, v_);
				const uint64_t lca_idx = phi(lca_idx_c);

				c = support_cbv.rank0(lca_idx_c) - support_cbv.rank0(b);
				return lca_idx;
			}

			// both pointers but not the same
			if(r_u != r_v) {
				uint64_t u_ = 0;
				if(r_u & 1) {
					u_ = r.select((r_u + 1) / 2 - 1);
				} else if(r_u == 0) u_ = u;
				else u_ = r.select(r_u / 2 - 1) + u - s.select(r_u - 1);

				uint64_t v_ = 0;
				if(r_v & 1) {
					v_ = r.select((r_v + 1) / 2 - 1);
				} else if(r_v == 0) v_ = v;
				else v_ = r.select(r_v / 2 - 1) + v - s.select(r_v - 1);

				if(v_ - u_ + 1 <= exc_sample_bs) {
					const int64_t exc_i = u_ - 2 * support_cbv.rank0(u_);
					auto min_exc = min_excess_scan(cbv, exc_i, u_, v_);
					c = support_cbv.rank0(min_exc.first) - support_cbv.rank0(b);
					const uint64_t lca_idx = phi(min_exc.first);
					return lca_idx;
				}
				auto [lca_idx_c, min_got] = min_excess(u_, v_);
				const uint64_t lca_idx = phi(lca_idx_c);

				c = support_cbv.rank0(lca_idx_c) - support_cbv.rank0(b);
				return lca_idx;
			}

			const uint64_t p = r_u;
			const uint64_t m = s.select(p - 1);
			const uint64_t rank_tree = ((p + 1) / 2) - 1;
			const uint64_t ptr_idx = this->p_s[rank_tree];
			const uint64_t ptr_idx_c = phi_1(ptr_idx);
			const uint64_t u_ = ptr_idx + u - m;
			const uint64_t v_ = ptr_idx + v - m;

#ifndef NDEBUG
			recursive_calls++;
#endif

			return lca(u_, v_, ptr_idx_c, c) - ptr_idx + m;
		}

		/**
		 * Given an index u inside the uncompressed balanced parenthesis sequence,
		 * return the inorder index of the node inside the tree. This function is 
		 * designed to be called after an lca(u) query.
		 * 
		 * @param a valid index u inside the original balanced parenthesis sequence
		 * @param the number of closing parentheses c between u and the starting position
		 * 		of the pointer where it lies (if it exists), 0 otherwise 
		 * @return the inorder index of the node u inside the tree
		 */
		inline uint64_t inorder(uint64_t u, uint64_t c = 0) {
			uint64_t r_u = s.rank(u);
			if(r_u & 1) {
				// open parenthesis of the pointer where u lies
				const uint64_t prev = s.select(r_u);

				// rank_0 up to that opening parenthesis
				const uint64_t prev_rank = rank0_outside(prev);

				assert(prev_rank + c >= 1);

				// minus 1 because indexes start from 0
				return prev_rank + c - 1;
			} else {
				return rank0_outside(u);
			}
		}

		scst& operator=(const scst& cbp) {
			cbv = pasta::BitVector(cbp.cbv.size(), 0);
			for(size_t i = 0; i < cbp.cbv.size(); i++)
				if(cbp.cbv[i])
					cbv[i] = 1;
			support_cbv = pasta::FlatRankSelect<pasta::OptimizedFor::ZERO_QUERIES>(cbv);
			innode_pointers_ones = cbp.innode_pointers_ones;
			std::vector< uint64_t > ones_innode_pointers;
			for(size_t i = 0; i < cbp.innode_pointers.rank(cbp.innode_pointers.size()); i++)
				ones_innode_pointers.push_back(cbp.innode_pointers.select(i));
			innode_pointers = sux::bits::EliasFano<>(ones_innode_pointers, cbp.innode_pointers.size());

			exc_samples = cbp.exc_samples;
			min_exc_samples = cbp.min_exc_samples;
			pos_min_exc_samples = cbp.pos_min_exc_samples;

			p_s = cbp.p_s;

			std::vector< uint64_t > ones_s;
			for(size_t i = 0; i < cbp.s.rank(cbp.s.size()); i++)
				ones_s.push_back(cbp.s.select(i));
			s = sux::bits::EliasFano<>(ones_s, cbp.s.size());

			std::vector< uint64_t > ones_r;
			for(size_t i = 0; i < cbp.r.rank(cbp.r.size()); i++)
				ones_r.push_back(cbp.r.select(i));
			r = sux::bits::EliasFano<>(ones_r, cbp.r.size());

			first_pointer_s = cbp.first_pointer_s;
			first_pointer_r = cbp.first_pointer_r;
			return *this;
		}

		inline uint64_t size_in_bits() const {
			return cbv.size() + precomputed_tables_size() + s.bitCount() + r.bitCount() +
				innode_pointers.bitCount() + 
				(sdsl::size_in_bytes(p_s) + support_cbv.space_usage() + 
				 sdsl::size_in_bytes(exc_samples) +
				 sdsl::size_in_bytes(min_exc_samples) + sdsl::size_in_bytes(pos_min_exc_samples) +
				 + sizeof(first_pointer_r) + 
				 sizeof(first_pointer_s)) * CHAR_BIT;
		}

		inline std::map<std::string, uint64_t> size_per_component() {
			std::map<std::string, uint64_t> sizes;
			sizes["01. pre comp tables"] = precomputed_tables_size();
			sizes["02. cbv"] = cbv.size();
			sizes["03. support_cbv"] = support_cbv.space_usage() * CHAR_BIT;
			sizes["04. innode pointers"] = (innode_pointers.bitCount()); 
			sizes["05. exc_samples"] = sdsl::size_in_bytes(exc_samples) * CHAR_BIT;
			sizes["06. min_exc_samples"] = sdsl::size_in_bytes(min_exc_samples) * CHAR_BIT;
			sizes["07. pos_min_exc_samples"] = sdsl::size_in_bytes(pos_min_exc_samples) * CHAR_BIT;
			sizes["08. s_support"] = (s.bitCount()); 
			sizes["09. r_support"] = (r.bitCount()); 
			sizes["10. pointers_s"] = sdsl::size_in_bytes(p_s) * CHAR_BIT;

			return sizes;
		}

		/**
		 * Given a position u of the uncompressed balanced parenthesis outside
		 * of a pointer, return the number of closed parentheses in the range [1,u),
		 * i.e., rank_{0}(B, u).
		 * 
		 * @param an index u of the original balanced parenthesis sequence outside of pointer
		 * @return the number of closed parentheses in the range [1,u), i.e. rank_{0}(B, u) 
		 */
		inline uint64_t rank0_outside(const uint64_t u) {

			const uint64_t v = phi_1(u);

			if(u > first_pointer_s) {
				const uint64_t rank = r.rank(v);
				const uint64_t len = u - first_pointer_s + 1;

				assert(v >= first_pointer_r);

				const uint64_t lenc = v - first_pointer_r + 1;

				assert(len > lenc);
				assert((len - lenc) % 2 == 0); 

				return (support_cbv.rank0(v) + (len - lenc)/2);
			} else {
				return support_cbv.rank0(v);
			}
		}

		inline uint64_t select0(const uint64_t i) {
			uint64_t l_in = 0;
			uint64_t r_in = innode_pointers_ones - 1;

			while(l_in < r_in) {
				uint64_t mid = (l_in + r_in + 1) / 2;
				if(innode_pointers.select(mid) < i) {
					l_in = mid;
				} else {
					r_in = mid - 1;
				}
			}

			const uint64_t search = i - innode_pointers.select(l_in);
			if(l_in + 1 == innode_pointers_ones) {
				const uint64_t pointer_s = s.select(l_in * 2 - 1);
				const uint64_t pointer_cbv = r.select(l_in - 1);

				const uint64_t inn_cbv = support_cbv.rank0(pointer_cbv + 1);
				const uint64_t inn_res = support_cbv.select0(search + inn_cbv);
				return pointer_s + inn_res - pointer_cbv;
			}

			uint64_t next_innode = innode_pointers.select(l_in + 1);
			uint64_t pointer_next = s.select((l_in + 1) * 2 - 1);
			uint64_t left_subtree = s.select((l_in + 1) * 2 - 2);
			if(i <= next_innode - ((pointer_next - left_subtree + 1) / 2)) {
				uint64_t pointer_cbv = 0;
				if(l_in > 0) {
					pointer_cbv = r.select(l_in - 1);
				}
				uint64_t inn_cbv = support_cbv.rank0(pointer_cbv + 1);
				uint64_t inn_res = support_cbv.select0(search + inn_cbv);
				return phi(inn_res);
			}

			// inside a pointer
			uint64_t acumm = next_innode - ((pointer_next - left_subtree + 1) / 2);
			uint64_t pos = left_subtree;
			while(i > acumm) {
				uint64_t new_l_in = s.rank(p_s[l_in]) / 2;
				uint64_t new_r_in = s.rank(p_s[l_in] + (pointer_next - left_subtree + 1)) / 2;
#ifndef NDEBUG
				recursive_calls++;
#endif

				const uint64_t pos_in_cbv = phi_1(p_s[l_in]);
				const uint64_t new_c =  support_cbv.rank0(pos_in_cbv) +
					(new_l_in == 0 ? 0 : innode_pointers.select(new_l_in) - support_cbv.rank0(r.select(new_l_in - 1) + 1));

				const uint64_t new_i = i - acumm;
				while(new_l_in < new_r_in) {
					uint64_t mid = (new_l_in + new_r_in + 1) / 2;
					if(innode_pointers.select(mid) < new_i + new_c) {
						new_l_in = mid;
					} else {
						new_r_in = mid - 1;
					}
				}

				const uint64_t curr_innode = innode_pointers.select(new_l_in);
				next_innode = innode_pointers.select(new_l_in + 1);
				pointer_next = s.select((new_l_in + 1) * 2 - 1);
				left_subtree = s.select((new_l_in + 1) * 2 - 2);
				// outside
				if(new_i + new_c <= next_innode - ((pointer_next - left_subtree + 1) / 2)) {

					const uint64_t pointer_cbv = (new_l_in == 0 ? 0 : r.select(new_l_in - 1));

					const uint64_t inn_cbv = support_cbv.rank0(pointer_cbv + 1);
					const uint64_t inn_res = support_cbv.select0(new_i + new_c - curr_innode + support_cbv.rank0(pointer_cbv + 1));
					return phi(inn_res) - p_s[l_in] + pos;
				}

				pos += left_subtree - p_s[l_in];
				l_in = new_l_in;
				acumm += next_innode - (pointer_next - left_subtree + 1) / 2 - new_c;
			}
			// this should never happen
			return -1;
		}

		void sample_excess_array(const std::vector<int64_t> &excess_array) {
			const size_t excess_size = excess_array.size();
			const size_t samples = (excess_size / exc_sample_bs) + 2;

			exc_samples = sdsl::int_vector<>(samples, 0);
			min_exc_samples = sdsl::int_vector<>(samples, 0);
			pos_min_exc_samples = sdsl::int_vector<>(samples, 0);

			for(size_t i = 0; i < excess_size; i += exc_sample_bs) {
				int64_t min_exc = excess_array[i];
				int64_t pos_min_exc = 0;

				for(size_t j = 1; i + j < excess_size && j < exc_sample_bs; ++j) {
					if(excess_array[i + j] < min_exc) {
						min_exc = excess_array[i + j];
						pos_min_exc = j;
					}
				}

				exc_samples[i / exc_sample_bs] = (i + exc_sample_bs - 1 < excess_size) ? excess_array[i + exc_sample_bs - 1] : excess_array.back();
				min_exc_samples[i / exc_sample_bs] = min_exc;
				pos_min_exc_samples[i / exc_sample_bs] = pos_min_exc;
			}

			sdsl::util::bit_compress(exc_samples);
			sdsl::util::bit_compress(min_exc_samples);
			sdsl::util::bit_compress(pos_min_exc_samples);
		}

		inline std::pair<size_t, int64_t> min_excess(const size_t l, const size_t h) const {

			if(h - l + 1 < exc_sample_bs) [[unlikely]] {
				//const auto exc_l = l - 2 * rank_cbv.rank0(l);
				const auto exc_l = l - 2 * support_cbv.rank0(l);
				return min_excess_scan(cbv, exc_l, l, h);
			}

			const size_t start = l / exc_sample_bs;
			const size_t end = h / exc_sample_bs;

			int64_t pos_min_exc = pos_min_exc_samples[start], min_exc = min_exc_samples[start];
			size_t sample_idx = start;

			for(size_t k = start + 1; k <= end; ++k) {
				if(min_exc_samples[k] < min_exc) {
					min_exc = min_exc_samples[k];
					pos_min_exc = pos_min_exc_samples[k];
					sample_idx = k;
				}
			}

			pos_min_exc = pos_min_exc + (sample_idx * exc_sample_bs);

			if(pos_min_exc >= l && pos_min_exc <= h) [[likely]] {
				return std::make_pair(pos_min_exc, min_exc);
			} else [[unlikely]] {

				//const auto exc_l = l - 2 * rank_cbv.rank0(l);
				const auto exc_l = l - 2 * support_cbv.rank0(l);

				// notice (start + 1) * exc_sample_bs - 1 is never >= cbv.size() because of the initial check
				std::tie(pos_min_exc, min_exc) = min_excess_scan(cbv, exc_l, l, (start + 1) * exc_sample_bs - 1);

				// todo: can be avoided
				for(size_t k = start + 1; k < end; ++k) {
					if(min_exc_samples[k] < min_exc) {
						min_exc = min_exc_samples[k];
						pos_min_exc = pos_min_exc_samples[k] + (k * exc_sample_bs);
						sample_idx = k;
					}
				}

				const auto [last_min_pos, last_min_exc] = min_excess_scan(cbv, exc_samples[end - 1], end * exc_sample_bs, h);
				pos_min_exc = (last_min_exc < min_exc) ? last_min_pos : pos_min_exc;
				min_exc = (last_min_exc < min_exc) ? last_min_exc : min_exc;

				return std::make_pair(pos_min_exc, min_exc);
			}
		}

		// information functions
		std::pair< double, double > average_subtree_size() {
			uint64_t subtree_sizes = 0;
			std::vector< uint64_t > sizes;
			for(size_t i = 0; i < p_s.size() * 2; i += 2) {
				sizes.push_back(s.select(i + 1) - s.select(i) + 1);
				subtree_sizes += sizes.back();
			}
			double average = (double) subtree_sizes / sizes.size();
			double des = 0;
			for(size_t i = 0; i < sizes.size(); i++) {
				des += (double) (sizes[i]  - average) * (sizes[i]  - average);
			}
			des /= (sizes.size() - 1);
			des = std::sqrt(des);
			return {average, des};
		}

		std::pair< double, double > average_length_chain()  {
			std::set< uint64_t > pointers_visited;
			std::vector< uint64_t > chains;

			for(size_t i = p_s.size(); i > 0; i--) {
				// already visited pointer
				if(pointers_visited.find(i) != pointers_visited.end()) continue;

				std::stack< std::pair< uint64_t, uint64_t >,
					std::vector< std::pair< uint64_t, uint64_t > > > next_pointer;

				next_pointer.push({i, 1});

				while(next_pointer.size() > 0) {
					auto [curr_i, curr_len] = next_pointer.top();
					next_pointer.pop();
					pointers_visited.insert(curr_i);

					uint64_t size_subtree = s.select(curr_i * 2  - 1) - s.select(curr_i * 2 - 2) + 1;
					uint64_t amount = (s.rank(p_s[curr_i - 1] + size_subtree) - s.rank(p_s[curr_i - 1])) / 2;
					uint64_t new_i = s.rank(p_s[curr_i - 1]) / 2;
					if(amount == 0) {
						chains.push_back(curr_len);
					} else {
						for(size_t j = 1; j <= amount; j++) {
							next_pointer.push({new_i + j, curr_len + 1});
						}
					}
				}
			}

			double average = 0;
			for(size_t i = 0; i < chains.size(); i++) {
				average += chains[i];
			}
			average /= chains.size();

			double des = 0;
			for(size_t i = 0; i < chains.size(); i++) {
				des += (chains[i] - average) * (chains[i] - average);
			}
			des /= (chains.size() - 1);
			des = std::sqrt(des);
			return {average, des};
		}
};
#endif

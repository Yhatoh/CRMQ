#pragma once

#include <cstdint>
#include <climits>
#include <vector>
#include <map>

#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>

#include "huffman.hpp"

template<typename T>
class huff_vector {

    using CanonicalHuffman = CanonicalHuffmanCode<T, std::vector<T>>;

    using SymbolTable = std::map<T, std::pair<uint64_t, uint64_t>>;

    using Code = std::pair<uint64_t, uint64_t>;

    CanonicalHuffman canonical_huffman;

    sdsl::bit_vector code_sequence;

    sdsl::sd_vector<> offsets;
    sdsl::select_support_sd<> select_offsets;

    uint64_t last_index;
    uint64_t sequence_length;

public:

    huff_vector() = default;

    explicit huff_vector(const std::vector<T> &data) : last_index(data.size() - 1) {

        std::map<T, uint64_t> frequencies;

        for(const T &d : data) frequencies[d]++;

        canonical_huffman = CanonicalHuffman(frequencies);

        SymbolTable symbols_to_codewords = canonical_huffman.enumerate_alphabet_code_pair(); 

        // compute the bit vector length

        sequence_length = 0;

        for(const auto &enc : symbols_to_codewords) {
            const T symbol = enc.first;
            const Code code = enc.second;
            const uint64_t occurrences = frequencies[symbol];
            sequence_length += code.first * occurrences;

            //std::cout << "symbol: " << symbol << " occ: " 
            //    << occurrences << " code: " << code.second << " len: " << code.first << std::endl;
        }

        //std::cout << "sequence len: " << sequence_length << std::endl;
        //std::cout << "last index: " << last_index << std::endl;

        code_sequence = sdsl::bit_vector(sequence_length, 0);

        std::vector<uint64_t> tmp_offsets;
        tmp_offsets.reserve(data.size());

        uint64_t offset = 0;

        for(const T &d : data) {
            tmp_offsets.push_back(offset);
            const Code code = symbols_to_codewords[d];
            const uint64_t len = code.first;
            for(uint64_t i = 0; i < len; ++i) { // notice: from lsb to msb
                code_sequence[offset + i] = (code.second & (1UL << i))? 1 : 0;
            }
            offset += len;
        }

        assert(tmp_offsets.size() == data.size());

        offsets = sdsl::sd_vector<>(tmp_offsets.begin(), tmp_offsets.end());
		sdsl::util::init_support(select_offsets, &offsets);

        /*for(int64_t i = 0; i < data.size(); ++i)
            std::cout << select_offsets(i + 1) << " ";
        std::cout << std::endl;*/
    }

    huff_vector& operator=(const huff_vector &hv) {
        canonical_huffman = hv.canonical_huffman;
        code_sequence = hv.code_sequence;
        offsets = hv.offsets;
        sdsl::util::init_support(select_offsets, &offsets);
        last_index = hv.last_index;
        sequence_length = hv.sequence_length;
		return *this;
	}

    const T operator[](uint64_t index) const {
        const uint64_t offset = select_offsets(index + 1);
        uint64_t next_offset;
        
        if(index == last_index) [[unlikely]]
            next_offset = sequence_length;
        else [[likely]]
            next_offset = select_offsets(index + 2);

        const uint64_t len = next_offset - offset;
        uint64_t code = 0UL;

        std::cout << "sequence len: " << sequence_length << std::endl;
        std::cout << "len: " << len << " start: " << offset << " end: " << next_offset << std::endl;

        // notice: from lsb to msb
        for(uint64_t i = offset, shift = 0UL; i < next_offset; ++i, ++shift) {
            code |= (code_sequence[i] << shift);
        }

        return canonical_huffman.decode(len, code);
    }

    inline uint64_t size() const {
        return code_sequence.size();
        /*return code_sequence.size() + (sdsl::size_in_bytes(offsets) + 
                + sdsl::size_in_bytes(select_offsets)) * CHAR_BIT;*/
        /*return canonical_huffman.size() + code_sequence.size() + (sdsl::size_in_bytes(offsets) + 
                + sdsl::size_in_bytes(select_offsets)) * CHAR_BIT;*/
    }
};
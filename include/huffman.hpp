#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <vector>
#include <algorithm>
#include <cassert>
#include <climits>

template <typename T>
class HuffmanTree {

    struct Node;

    using NodePtr = std::shared_ptr<Node>;

    struct Node {
        T key;
        uint64_t freq;
        NodePtr next, left, right;
        Node(uint64_t freq) : key(), freq(freq), next(), left(), right() {}
        Node(T key, uint64_t freq)
            : key(key), freq(freq), next(), left(), right() {}
    };

    uint64_t number_of_alphabet;
    NodePtr root;

public:

    HuffmanTree() = default;

    explicit HuffmanTree(const std::map<T, uint64_t> &alphabet_count) : number_of_alphabet(alphabet_count.size()) {
        std::vector<NodePtr> leaves;
        leaves.reserve(number_of_alphabet);
        for (auto &&[key, freq] : alphabet_count) {
            leaves.push_back(std::make_shared<Node>(key, freq));
        }
        std::sort(leaves.begin(), leaves.end(),
                [](const NodePtr &a, const NodePtr &b) -> bool {
                    return a->freq < b->freq;
                });
        for (int i = 0; i < number_of_alphabet - 1; i++) {
            leaves[i]->next = leaves[i + 1];
        }
        NodePtr insert_pos = leaves[0];
        NodePtr look_pos = leaves[0];
        while (look_pos->next != nullptr) {
            NodePtr n =
                std::make_shared<Node>(look_pos->freq + look_pos->next->freq);
            n->left = look_pos;
            n->right = look_pos->next;
            while (insert_pos->next != nullptr &&
                insert_pos->next->freq <= n->freq) {
                insert_pos = insert_pos->next;
            }
            n->next = insert_pos->next;
            insert_pos->next = n;
            look_pos = look_pos->next->next;
        }
        root = look_pos;

        auto delete_next = [](auto self, NodePtr n) -> void {
            n->next = nullptr;
            if (n->left != nullptr) {
                self(self, n->left);
            }
            if (n->right != nullptr) {
                self(self, n->right);
            }
        };
        delete_next(delete_next, root);
    }

    std::vector<std::pair<T, uint64_t>> compute_length() {
        std::vector<std::pair<T, uint64_t>> res;
        res.reserve(number_of_alphabet);

        auto dfs = [&](auto self, NodePtr n, uint64_t d = 0) -> void {
            if (n->left == nullptr) {
                res.push_back({n->key, d});
            } else {
                self(self, n->left, d + 1);
                self(self, n->right, d + 1);
            }
        };
        dfs(dfs, root);

        return res;
    }
};

template<typename T, typename AlphabetArray>
class CanonicalHuffmanCode {

    AlphabetArray alphabets;
    std::vector<uint64_t> first_index, first_code;
    uint64_t number_of_alphabet, maximum_code_length;

public:

    CanonicalHuffmanCode() = default;

    explicit CanonicalHuffmanCode(const std::map<T, uint64_t> &alphabet_count) : number_of_alphabet(alphabet_count.size()) {
        HuffmanTree<T> ht(alphabet_count);
        auto alphabet_and_code_lengths = ht.compute_length();
        sort(alphabet_and_code_lengths.begin(), alphabet_and_code_lengths.end(),
            [](const std::pair<T, uint64_t> &a, const std::pair<T, uint64_t> &b)
                -> bool { return a.second < b.second; });

        std::vector<T> alphabets_vec;
        alphabets_vec.reserve(number_of_alphabet);
        for (auto &&[alphabet, code_length] : alphabet_and_code_lengths) {
            alphabets_vec.push_back(alphabet);
        }
        alphabets = std::move(alphabets_vec);

        maximum_code_length = alphabet_and_code_lengths.back().second;
        assert(maximum_code_length <= 64);

        std::vector<uint64_t> codes;
        codes.reserve(number_of_alphabet);
        codes.push_back(0);
        for (int i = 1; i < number_of_alphabet; i++) {
            codes.push_back((codes.back() + 1)
                            << (alphabet_and_code_lengths[i].second -
                                alphabet_and_code_lengths[i - 1].second));
        }

        first_index.resize(maximum_code_length + 1, number_of_alphabet);
        first_code.resize(maximum_code_length + 1);
        for (int i = number_of_alphabet - 1; i >= 0; i--) {
            first_index[alphabet_and_code_lengths[i].second] = i;
            first_code[alphabet_and_code_lengths[i].second] = codes[i];
        }
        for (int l = maximum_code_length - 1; l >= 1; l--) {
            if (first_index[l] == number_of_alphabet) {
                first_index[l] = first_index[l + 1];
                first_code[l] = first_code[l + 1] >> 1;
            }
        }
    }

    CanonicalHuffmanCode& operator=(const CanonicalHuffmanCode &canonical_huff) {
        alphabets = canonical_huff.alphabets;
        first_index = canonical_huff.first_index;
        first_code = canonical_huff.first_code;
        number_of_alphabet = canonical_huff.number_of_alphabet;
        maximum_code_length = canonical_huff.maximum_code_length;
        return *this;
    }

    std::map<T, std::pair<uint64_t, uint64_t>> enumerate_alphabet_code_pair() const {
        std::map<T, std::pair<uint64_t, uint64_t>> result;
        uint64_t code_length = 0, code = 0;
        for (int i = 0; i < number_of_alphabet; i++) {
            while (code_length < maximum_code_length &&
                first_code[code_length + 1] <= (code << 1)) {
                code_length++;
                code <<= 1;
            }
            result[alphabets[i]] = {code_length, code};
            code++;
        }
        return result;
    }

    T decode(uint64_t code_length, uint64_t code) const {
        return alphabets[first_index[code_length] + code - first_code[code_length]];
    }

    uint64_t size() const {
        return ((first_index.size() + first_code.size()) * sizeof(uint64_t) +
                + alphabets.size() * sizeof(alphabets)) * CHAR_BIT;
    }
};
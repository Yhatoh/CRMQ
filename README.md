# 🗜️ Compressed Range Minimum Queries

This repository implements a compressed encoding for the Range Minimum Query (RMQ) problem. Our method achieves compression by exploiting maximal identical subtrees in the Cartesian tree of the underlying input array, obtaining its minimal Directed Acyclic Graph (DAG) representation. The construction is performed in linear time and space, directly from a succinct representation of the tree.

Experimental evaluation on 11 real-world datasets shows that the proposed approach is highly space-efficient, requiring as little as $0.11n$ bits in practice, while delivering fast query times.

## Building the project

Clone the repository.

```bash
git clone --recursive git@github.com:Yhatoh/CRMQ.git
cd cct
```

Build the `libsais`, `bit_vector`, and `sdsl-lite` libraries separately, for example:

```bash
cd libsais
cmake .
make -j8
```

Then build the `cct` library.

```bash
mkdir build
cd build
cmake ..
make -j8
```

## Reproducing our results

You can find our benchmarking suite and instructions on how to reproduce our results
at the following links:

* For uncompressed RMQ solutions: https://github.com/FilippoLari/RMQ-experiments
* For compressed RMQ solutions: https://github.com/Yhatoh/CRMQ-Experiments

## Credits

If you use this project in a scientific article, please cite the following paper:

```
@inproceedings{CarmonaL26,
    author = {Gabriel Carmona and Filippo Lari},
    title = {Compressing Highly Repetitive Binary Trees with an Application to Range Minimum Queries},
    booktitle = {Proc. 24th Symposium on Experimental Algorithms},
    year = {2026},
}
```

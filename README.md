# Compressed Range Minimum Queries

This repository implements the Subtree Compressed Succinct Tree and they use solving
range mininum queries.

## Building the project

Clone the repository.

```bash
git clone --recursive git@github.com:Yhatoh/CRMQ.git
cd cct
```

Build the `libsais`, `bit_vector` and `sdsl-lite` libraries separately, e.g.

```bash
cd libsais
cmake .
make -j8
```

Then build the `cct` library.

```
mkdir build
cd build
cmake ..
make -j8
```

## Reproducing our results

You can find our benchmarking suit and instructions on how to reproduce our results
at the following links:

* For non-compressed RMQs solutions: [NON COMPRESSED RMQ-EXPERIMENTS](https://github.com/FilippoLari/RMQ-experiments)
* For compressed RMQs solutions: [COMPRESSED RMQ-EXPERIMENTS](https://github.com/Yhatoh/CRMQ-Experiments)

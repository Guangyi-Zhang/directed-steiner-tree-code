
# Directed Steiner Tree

Project structure

- `dst`: code for the algorithms
- `main`: code for experimental evaluation
- `test`: code for test cases
- `example`: a simple example for demonstration

The project requires the following prerequisites

- C++ 17 or later
- Boost
- CMake

In Ubuntu, these can be installed via

```bash
apt install -y make \
               build-essential \
               cmake \
               libboost-all-dev
```


## Example

See `./example` for a concrete example.
```bash
# or, bash example.sh
cmake -S ./ -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_TEST=0 -DUSE_MAIN=0 -DUSE_EXAMPLE=1
cmake --build build
./build/example/Example
```

For input as follows,
```cpp
std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                       std::make_pair(0,2), 
                                       std::make_pair(0,3), 
                                       std::make_pair(0,4), 
                                       std::make_pair(1,11), 
                                       std::make_pair(2,12), 
                                       std::make_pair(3,13), 
                                       std::make_pair(4,11), 
                                       std::make_pair(4,12), 
                                       std::make_pair(4,13)};
std::vector<double> weights {1,1,1,2.5, 2,2,2, 1,1,1};
std::vector<int> terms {11,12,13};
int root = 0;
```

The output with the `fast_level2` algorithm is 
```
cost 5.5, with time 0.014 ms
└──0
    └──4 (002.5)
        ├──13 (001.0)
        │   └──16 (000.0)
        ├──12 (001.0)
        │   └──15 (000.0)
        └──11 (001.0)
            └──14 (000.0)
```

We refer the user to `./test` and `./main` for more usages and algorithms.

## Tests

```bash
cmake -S ./ -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_TEST=1 -DUSE_MAIN=0 -DUSE_EXAMPLE=0
cmake --build build
./build/test/Test
```

## Datasets

- https://networkrepository.com/advogato.php
- https://snap.stanford.edu/data/soc-Epinions1.txt.gz
- https://snap.stanford.edu/data/web-Google.txt.gz
- https://snap.stanford.edu/data/wiki-topcats.html
- https://snap.stanford.edu/data/soc-LiveJournal1.html

<!-- 
- https://snap.stanford.edu/data/soc-pokec-relationships.txt.gz

- https://snap.stanford.edu/data/cit-HepPh.txt.gz
- https://lfs.aminer.cn/lab-datasets/citation/DBLP_citation_Sep_2013.rar

- https://snap.stanford.edu/data/roadNet-CA.txt.gz
- https://snap.stanford.edu/data/ERC20-stablecoins.zip

- https://users.cs.utah.edu/~lifeifei/research/tpq/SF.cedge
- https://users.cs.utah.edu/~lifeifei/research/tpq/cal.cedge
-->

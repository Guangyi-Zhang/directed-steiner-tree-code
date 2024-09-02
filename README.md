
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

In Ubuntu, these can be installed by

```bash
apt install -y make \
               build-essential \
               cmake \
               libboost-all-dev
```


## Example

See `./example` for a concrete example.
Users are also referred to `./test` and `./main` for more usages.

```bash
# or, base example.sh
cmake -S ./ -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_TEST=0 -DUSE_MAIN=0 -DUSE_EXAMPLE=1
cmake --build build
./build/example/Example
```

Output:
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

## Tests

```bash
cmake -S ./ -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_TEST=1 -DUSE_MAIN=0
cmake --build build
./build/test/Test
```

## Datasets

- https://networkrepository.com/advogato.php
- https://snap.stanford.edu/data/soc-Epinions1.txt.gz
- https://snap.stanford.edu/data/web-Google.txt.gz
- https://snap.stanford.edu/data/soc-pokec-relationships.txt.gz
- https://snap.stanford.edu/data/soc-LiveJournal1.html

<!-- 
- https://snap.stanford.edu/data/cit-HepPh.txt.gz
- https://lfs.aminer.cn/lab-datasets/citation/DBLP_citation_Sep_2013.rar

- https://snap.stanford.edu/data/roadNet-CA.txt.gz
- https://snap.stanford.edu/data/ERC20-stablecoins.zip

- https://users.cs.utah.edu/~lifeifei/research/tpq/SF.cedge
- https://users.cs.utah.edu/~lifeifei/research/tpq/cal.cedge
-->
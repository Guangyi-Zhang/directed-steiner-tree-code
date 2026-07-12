
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
               cmake
```

and locally under `third_party/` (auto-detected by CMake):

```bash
apt download libboost1.83-dev
dpkg-deb -x libboost1.83-dev_*.deb ./third_party/boost-dev
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
The main experiment binary accepts `level1`, `adaptive_level1`, `level2`,
`adaptive_level2` (also `root_adaptive_level2`), `adaptive_fast_level2`
(also `root_adaptive_fast_level2`), `full_adaptive_level2`,
`full_adaptive_fast_level2`, `fast_level2`, `level3`, and `fast_level3`.
`adaptive_fast_level2` uses the dynamic-prefix root-adaptive implementation:
terminal-side distances are discovered progressively and reused while each
greedy iteration refreshes only the root-side distances.
For exact adaptive 2-level methods, each greedy iteration keeps the top 5
scored-density candidates by default and reranks that shortlist by realized
marginal density.  Fast adaptive 2-level methods certify only the top scored
candidate, then rerank a bounded pool of concrete observed prefixes discovered
during that certificate.  Scored pools are ordered by lower scored density,
fewer terminals, then lower center id; realized reranking uses lower realized
density with the same scored tie keys. Density ordering uses exact raw values;
`EPSILON` comparisons are reserved for assertions and diagnostics. Pass
`rerank_L = 1` to the C++ API to recover single-candidate scored selection.
All adaptive 2-level C++ entry points reject `rerank_L < 1` with
`std::invalid_argument`. Another terminal is added only when it strictly
improves scored density.

Fast adaptive 2-level results also report `observed_pool_consider_calls`,
`observed_pool_key_changes`, `observed_pool_admissions`,
`observed_pool_rejections`, `observed_pool_materializations`, and
`observed_pool_invariant_violations`. These 64-bit counters distinguish cheap
snapshot maintenance from final `PartialTree` materialization. Normally,
materializations are at most `rerank_L` per greedy iteration and invariant
violations are zero.


## Tests

```bash
cmake -S ./ -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_TEST=1 -DUSE_MAIN=0 -DUSE_EXAMPLE=0
cmake --build build
./build/test/Test --source-file-exclude='*test_level3.cpp'
```

## Datasets

- https://networkrepository.com/advogato.php
- https://snap.stanford.edu/data/soc-Epinions1.txt.gz
- https://snap.stanford.edu/data/web-Google.txt.gz
- https://snap.stanford.edu/data/wiki-topcats.html
- https://snap.stanford.edu/data/soc-LiveJournal1.html
- https://mcauleylab.ucsd.edu/public_datasets/data/amazon_v2/categoryFilesSmall/all_csv_files.csv
- http://konect.cc/networks/dbpedia-link/

<!-- 
- https://snap.stanford.edu/data/soc-pokec-relationships.txt.gz

- https://snap.stanford.edu/data/cit-HepPh.txt.gz
- https://lfs.aminer.cn/lab-datasets/citation/DBLP_citation_Sep_2013.rar

- https://snap.stanford.edu/data/roadNet-CA.txt.gz
- https://snap.stanford.edu/data/ERC20-stablecoins.zip

- https://users.cs.utah.edu/~lifeifei/research/tpq/SF.cedge
- https://users.cs.utah.edu/~lifeifei/research/tpq/cal.cedge
-->

For synthetic datasets, see `data.py`, and run `uv sync` to install the dependencies.

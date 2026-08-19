
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

The `example/` target demonstrates the C++ API on a small directed graph.
See `example/src/main.cpp` for the full program.

### Build and run

```bash
cmake -S ./ -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_TEST=0 -DUSE_MAIN=0 -DUSE_EXAMPLE=1
cmake --build build
./build/example/Example fast_level2
```

Or use the helper script (pass the method name as the first argument):

```bash
bash example.sh fast_level2
```

The example binary accepts `level2` and `fast_level2`.
The experiment binary in `main/` exposes the full method list via `--method`.

### Code sketch

```cpp
#include "dst/dst.hpp"

// Input graph: root = 0, terminals = {11, 12, 13}
//   0→1, 0→2, 0→3  (weight 1)
//   0→4             (weight 2.5)
//   1→11, 2→12, 3→13  (weight 2)
//   4→11, 4→12, 4→13  (weight 1)
std::vector<std::pair<int,int>> edges = {{0,1}, {0,2}, {0,3}, {0,4},
                                         {1,11}, {2,12}, {3,13},
                                         {4,11}, {4,12}, {4,13}};
std::vector<float> weights = {1, 1, 1, 2.5, 2, 2, 2, 1, 1, 1};
std::vector<int> terms = {11, 12, 13};
int root = 0;

dst::DST dt(edges, weights, root, terms);
auto partree = dt.fast_level2_alg();
auto tree = partree->to_tree();
```

The `DST` constructor removes unreachable vertices, then attaches zero-weight
dummy edges so each terminal can be tracked internally. Every returned tree is
checked for feasibility (a root-directed path to each terminal) before printing.

### Output

Running `fast_level2` on this instance gives cost **5.5**:

```
cost 5.5, with time 0.023 ms
└──0
    └──4 (002.5)
        ├──13 (001.0)
        │   └──16 (000.0)
        ├──12 (001.0)
        │   └──15 (000.0)
        └──11 (001.0)
            └──14 (000.0)
```

Edge weights are shown in parentheses. Dummy terminal nodes used internally are also printed with zero weight.

### Algorithms

**Level 1** — direct root-to-terminal paths.

| Method | C++ API | Summary |
|--------|---------|---------|
| `level1` | `level1_alg()` | One multi-terminal Dijkstra from the root; take the union of shortest paths to all terminals. |
| `adaptive_level1` | `adaptive_level1_alg()` | Greedy variant: repeatedly connect the nearest uncovered terminal, then zero out covered edges before the next Dijkstra. |

**Level 2** — root → center → terminal bunches, chosen greedily by terminal density.

| Method | C++ API | Summary |
|--------|---------|---------|
| `level2` | `level2_alg()` | Exact 2-level greedy: each step materializes every center candidate and picks the best realized density. |
| `fast_level2` | `fast_level2_alg()` | Same greedy rule as `level2`, but prunes centers with coordinated reverse Dijkstra and cheap lower bounds. |
| `adaptive_level2` | `adaptive_level2_alg(rerank_L=5)` | Root-adaptive 2-level: terminal-side distances are fixed; root-side distances are refreshed each iteration. Keeps the top `rerank_L` scored candidates and reranks them by realized density. Alias: `root_adaptive_level2`. |
| `full_adaptive_level2` | `full_adaptive_level2_alg(rerank_L=5)` | Fully adaptive 2-level: both root-side and terminal-side distances are refreshed each iteration, with the same top-`rerank_L` reranking. |
| `adaptive_fast_level2` | `adaptive_fast_level2_alg(rerank_L=5)` | Fast root-adaptive 2-level with dynamic-prefix center indices (Fenwick trees). Alias: `root_adaptive_fast_level2`. |
| `adaptive_fast_level2_no_fenwick` | `adaptive_fast_level2_no_fenwick_alg(rerank_L=5)` | Same search as `adaptive_fast_level2`, but uses scan-based center indices (less auxiliary memory, linear prefix queries). Alias: `root_adaptive_fast_level2_no_fenwick`. |
| `full_adaptive_fast_level2` | `full_adaptive_fast_level2_alg(rerank_L=5)` | Fully adaptive 2-level with fast pruning and observed-prefix reranking. |

Adaptive 2-level entry points require `rerank_L >= 1` and throw `std::invalid_argument` otherwise. Pass `rerank_L = 1` to disable multi-candidate reranking.


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

#pragma once

#include <vector>
#include <limits>

namespace dst {

  // Compressed Sparse Row graph with both forward and reverse adjacency lists.
  // fwd_* stores outgoing edges; rev_* stores incoming edges for reverse traversal.
  // Offset arrays have size n+1, where offset[i] to offset[i+1] marks the range
  // of edges for node i in the target/weight arrays.
  struct CSRGraph {
    int n = 0;  // number of nodes
    int m = 0;  // number of edges
    std::vector<int> fwd_offset, fwd_target;
    std::vector<float> fwd_weight;
    std::vector<int> rev_offset, rev_target;
    std::vector<float> rev_weight;

    // Returns weight of edge (u,v), or max float if edge doesn't exist.
    float edge_weight(int u, int v) const {
      for (int i = fwd_offset[u]; i < fwd_offset[u+1]; i++)
        if (fwd_target[i] == v) return fwd_weight[i];
      return std::numeric_limits<float>::max();
    }

    // Sets edge weight to 0 in both forward and reverse directions.
    void zero_edge(int u, int v) {
      for (int i = fwd_offset[u]; i < fwd_offset[u+1]; i++) {
        if (fwd_target[i] == v) { fwd_weight[i] = 0.0f; break; }
      }
      for (int i = rev_offset[v]; i < rev_offset[v+1]; i++) {
        if (rev_target[i] == u) { rev_weight[i] = 0.0f; break; }
      }
    }
  };

}  // namespace dst

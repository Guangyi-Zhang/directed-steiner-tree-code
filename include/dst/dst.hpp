#pragma once

#include <string>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <algorithm>


namespace dst {

  constexpr int NONVERTEX {-1};
  
  std::pair<std::unordered_map<int,double>, 
            std::unordered_map<int,int>> dijkstra(std::map<int, std::vector<int>> adj, 
                                                  std::map<std::pair<int,int>, double> edgeweight, 
                                                  int source,
                                                  bool reverse=false) {
    std::unordered_map<int,double> distances;
    std::unordered_map<int,int> trace;
    std::priority_queue<std::tuple<double,int,int>, 
                        std::vector<std::tuple<double,int,int>>, 
                        std::greater<std::tuple<double,int,int>>> pq;
    distances[source] = 0;
    pq.emplace(0, NONVERTEX, source);

    while (!pq.empty()) {
      double du;
      int u_prev, u; 
      std::tie(du, u_prev, u) = pq.top();
      pq.pop();

      if (distances.find(u) != distances.end() and du > distances[u]) 
        continue;

      trace[u] = u_prev;
      for (const auto& v: adj[u]) {
          double weight = reverse? edgeweight[{v,u}] : edgeweight[{u,v}];

          if (distances.find(v) == distances.end() or distances[v] > distances[u] + weight) {
              distances[v] = distances[u] + weight;
              pq.emplace(distances[v], u, v);
          }
      }
    }

    return std::make_pair(distances, trace);
  }


  /**
   * @brief A class for saying hello in multiple languages
   */
  class DST {
  public:
    std::string name;

    int root;
    std::vector<int> terms;
    std::vector<int> terms_dm; // dummy
    std::map<int, int> terms_map;
    std::map<std::pair<int,int>, double> w;
    std::map<int, std::vector<int>> adj;
    std::map<int, std::vector<int>> adj_r; // reverse adj
    std::unordered_set<int> V; // excluding dummy terminals

    DST(std::vector<std::pair<int,int>> edges, std::vector<double> edgeweights, int root, std::vector<int> terms) :
        root {root},
        terms {terms}
    {
        int v_max = 0;
        for (size_t i = 0; i < edges.size(); i++) {
            auto p = edges[i];
            adj[p.first].push_back(p.second);
            adj_r[p.second].push_back(p.first);
            w[{p.first, p.second}] = edgeweights[i];
            V.insert(p.first);
            V.insert(p.second);
            v_max = std::max({p.first, p.second, v_max});
        }

        // create dummy terminals
        int i {1};
        for (auto t: terms) {
          int t_dm = v_max + i++;
          terms_dm.push_back(t_dm);
          terms_map[t] = t_dm;
          terms_map[t_dm] = t;

          adj[t].push_back(t_dm);
          adj_r[t_dm].push_back(t);
          w[{t, t_dm}] = 0;
        }
    }

    double naive_alg() {
      auto p_ = dijkstra(adj, w, root);
      auto dists = p_.first;
      auto trace = p_.second;
      std::set<std::pair<int,int>> edges_marked;
      for (auto t: terms_dm) {
        if (trace.find(t) == trace.end())
            continue; // disconnected graph

        int v {t};
        while (trace[v] != NONVERTEX) {
          edges_marked.insert({trace[v], v});
          v = trace[v];
        }
      }

      double total {0};
      for (auto e: edges_marked)
        total = total + w[{e.first, e.second}];
      return total;
    }

  };

}  // namespace greeter

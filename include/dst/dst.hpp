#pragma once

#include <limits>
#include <utility>
#include <numeric> // std::iota
#include <algorithm>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <tuple>


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
      double d_u;
      int u_prev, u; 
      std::tie(d_u, u_prev, u) = pq.top();
      pq.pop();

      if (distances.find(u) != distances.end() and d_u > distances[u]) 
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


  template <typename T>
  std::vector<size_t> argsort(const std::vector<T> &v) {
    // https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0); // fill idx by increasing values

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    stable_sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
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


    std::pair<double, std::vector<int>> 
        level2_rooted_at_v(int v, 
                           double d_rv, 
                           std::unordered_set<int> terms_left,
                           std::unordered_map<int, std::unordered_map<int,double>> dists_t) {
      // collect all t's for the current v
      std::vector<int> terms_left_v;
      std::vector<double> ds_vt;
      for (auto t: terms_left) {
        if (dists_t[t].find(v) == dists_t[t].end()) 
          continue; // t is not reachable from v
        
        terms_left_v.push_back(t);
        ds_vt.push_back(dists_t[t][v]);
      }

      // sort distances from v to t's
      double den_v = std::numeric_limits<double>::max();
      std::vector<int> terms_cov_v;
      auto idxs = argsort(ds_vt);
      for (size_t i=0; i < ds_vt.size(); i++) {
        auto idx = idxs[i];
        double den_i = (d_rv + ds_vt[idx]) / (i+1);
        if (den_i < den_v) {
          den_v = den_i;
          terms_cov_v.push_back(terms_left_v[idx]);
        }
        else {
          break; // stop once adding a t increases den_v
        }
      }

      return std::make_pair(den_v, terms_cov_v);
    }


    std::pair<double, std::unordered_map<int,int>> level2_alg() {
      // dijktra from the root
      auto p_r = dijkstra(adj, w, root);
      auto dists_r = p_r.first;
      auto trace_r = p_r.second;

      // dijktra from each terminal
      std::unordered_map<int, std::unordered_map<int,int>> trace_t;
      std::unordered_map<int, std::unordered_map<int,double>> dists_t;
      for (auto t: terms_dm) {
        auto p_t = dijkstra(adj_r, w, t, true);
        dists_t[t] = p_t.first;
        trace_t[t] = p_t.second;
      }

      // iteratively add 2-level partial trees
      std::unordered_map<int,int> trace; // trace "short-cut" chosen 2-level trees
      std::set<std::pair<int,int>> edges_marked;
      std::unordered_set<int> terms_left(terms_dm.begin(), terms_dm.end());
      while (terms_left.size() > 0) {
        // used to keep track of the best v so far
        int v_min {NONVERTEX};
        std::vector<int> terms_min;
        double den_min = std::numeric_limits<double>::max();

        // enum all v as the middle vertex in a 2-level tree
        for (auto v: V) {
          double d_rv {dists_r[v]};
          auto p_v = level2_rooted_at_v(v, d_rv, terms_left, dists_t);
          auto den_v = p_v.first;
          auto terms_cov_v = p_v.second;

          // keep the best across all v
          if (den_v < den_min) {
            v_min = v;
            terms_min = terms_cov_v;
            den_min = den_v;
          }
        }

        // the rest terminals are not reachable
        if (terms_min.size() == 0)
          break;

        // merge the best 2-level partial tree
        if (trace.find(v_min) == trace.end())
          trace[v_min] = root;
        for (auto t: terms_min) {
          trace[t] = v_min;
          terms_left.erase(t);
        }
        // mark the path from root to v
        auto u = v_min;
        while (trace_r[u] != NONVERTEX) {
          edges_marked.insert({trace_r[u], u});
          u = trace_r[u];
        }
        // mark the path from v to each t
        for (auto t: terms_min) {
          u = v_min;
          while (trace_t[t][u] != NONVERTEX) {
            edges_marked.insert({u, trace_t[t][u]});
            u = trace_t[t][u];
          }
        }
      }

      double total {0};
      for (auto e: edges_marked) {
        total += w[{e.first, e.second}];
      }
      return std::make_pair(total, trace);
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

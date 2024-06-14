#pragma once

#include <limits>
#include <type_traits> // is_same
#include <utility>
#include <cassert>
#include <numeric> // std::iota
#include <algorithm>
#include <iterator> // inserter
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <tuple>

#include <fmt/ranges.h>


namespace dst {

  constexpr int NONVERTEX {-1};


  template <typename T, typename U>
  bool has_key(const T& container, const U& key) {
    // only works for set and map!
    return (container.find(key) != container.end());
  }

  
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

      if (has_key(distances, u) and d_u > distances[u]) 
        continue;

      trace[u] = u_prev;
      for (const auto& v: adj[u]) {
          double weight = reverse? edgeweight[{v,u}] : edgeweight[{u,v}];

          if (not has_key(distances, v) or distances[v] > distances[u] + weight) {
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


  class PartialTree {
  public:
    double cost, cost_sc;
    std::unordered_map<int,int> trace_sc; // trace "short-cut" chosen 2-level trees
    std::unordered_map<int,int> trace;
    std::unordered_set<int> terms_cov;

    PartialTree(int root) {
      trace[root] = NONVERTEX;
      trace_sc[root] = NONVERTEX;
    }

    double density() const {
      return cost / terms_cov.size();
    }

    bool zero_coverage() const {
      return terms_cov.size() == 0 ? true: false;
    }

    std::set<std::pair<int,int>> edges() const {
      std::set<std::pair<int,int>> es;

      for (auto t: terms_cov) {
        auto v {t};
        while (trace.at(v) != NONVERTEX) {
          es.insert(std::make_pair(trace.at(v), v));
          v = trace.at(v);
        }
      }

      return es;
    }
  };


  /**
   * @brief A class for saying hello in multiple languages
   */
  class DST {
  public:
    std::string name;

    int root;
    std::vector<int> terms;
    std::unordered_set<int> terms_dm;// dummy
    std::map<int, int> terms_map;
    std::map<std::pair<int,int>, double> w;
    std::map<int, std::vector<int>> adj;
    std::map<int, std::vector<int>> adj_r; // reverse adj
    std::unordered_set<int> V; // excluding dummy terminals

    DST(  std::vector<std::pair<int,int>> edges, 
          std::vector<double> edgeweights, 
          int root, 
          std::vector<int> terms) :
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
          terms_dm.insert(t_dm);
          terms_map[t] = t_dm;
          terms_map[t_dm] = t;

          adj[t].push_back(t_dm);
          adj_r[t_dm].push_back(t);
          w[{t, t_dm}] = 0;
        }
    }


    std::pair<double, std::vector<int>> 
        level2_through_v(
          int v, 
          double d_rv, 
          std::unordered_map<int, std::unordered_map<int,double>> dists_t,
          std::unordered_set<int> terms_left
        ){
      // collect all t's for the current v
      std::vector<int> terms_left_v;
      std::vector<double> ds_vt;
      for (auto t: terms_left) {
        if (not has_key(dists_t[t], v)) 
          continue; // t is not reachable from v
        
        terms_left_v.push_back(t);
        ds_vt.push_back(dists_t[t][v]);
      }

      // sort distances from v to t's
      double cost_v {d_rv};
      std::vector<int> terms_cov_v;
      auto idxs = argsort(ds_vt);
      for (size_t i=0; i < ds_vt.size(); i++) {
        auto idx = idxs[i];
        double den_i = (cost_v + ds_vt[idx]) / (i+1);
        if (i == 0 or den_i < cost_v / i) {
          cost_v += ds_vt[idx];
          terms_cov_v.push_back(terms_left_v[idx]);
        }
        else {
          break; // stop once adding a t increases den_v
        }
      }

      return std::make_pair(cost_v, terms_cov_v);
    }


    PartialTree level2_rooted_at_r(
        int r,
        std::unordered_set<int> V_cand,
        std::unordered_set<int> terms_cand
    ) {
      // dijktra from the root
      auto p_r = dijkstra(adj, w, r);
      auto dists_r = p_r.first;
      auto trace_r = p_r.second;

      // dijktra from each terminal
      std::unordered_map<int, std::unordered_map<int,int>> trace_t;
      std::unordered_map<int, std::unordered_map<int,double>> dists_t;
      for (auto t: terms_cand) {
        auto p_t = dijkstra(adj_r, w, t, true);
        dists_t[t] = p_t.first;
        trace_t[t] = p_t.second;
      }

      // iteratively add 2-level partial trees
      PartialTree par {root};
      std::unordered_set<int> terms_left(terms_cand.begin(), terms_cand.end());
      while (terms_left.size() > 0) {
        // used to keep track of the best v so far
        int v_min {NONVERTEX};
        std::vector<int> terms_min;
        double den_min = std::numeric_limits<double>::max();

        // enum all v as the middle vertex in a 2-level tree
        for (auto v: V_cand) {
          double d_rv {dists_r[v]};
          auto p_v = level2_through_v(v, d_rv, dists_t, terms_left);
          auto cost_v = p_v.first;
          auto terms_cov_v = p_v.second;

          // keep the best across all v
          if (cost_v/terms_cov_v.size() < den_min) {
            v_min = v;
            terms_min = terms_cov_v;
            den_min = cost_v/terms_cov_v.size();
          }
        }

        // the rest terminals are not reachable
        if (terms_min.size() == 0)
          break;

        // merge the best 2-level partial tree
        if (not has_key(par.trace_sc, v_min)) {
          par.trace_sc[v_min] = root;
          par.cost_sc += dists_r[v_min];
        }
        for (auto t: terms_min) {
          par.trace_sc[t] = v_min;
          par.cost_sc += dists_t[t][v_min];
          terms_left.erase(t);
        }
        // mark the path from root to v
        auto u = v_min;
        while (trace_r[u] != NONVERTEX) {
          par.trace[u] = trace_r[u];
          u = trace_r[u];
        }
        // mark the path from v to each t
        for (auto t: terms_min) {
          u = v_min;
          while (trace_t[t][u] != NONVERTEX) {
            par.trace[trace_t[t][u]] = u;
            u = trace_t[t][u];
          }
        }
      }

      for (auto t: terms_cand) {
        if (not has_key(terms_left, t))
          par.terms_cov.insert(t);
      }
      double total {0};
      for (auto e: par.edges()) {
        total += w[{e.first, e.second}];
      }
      par.cost = total;
      return par;
    }


    PartialTree level2_alg() {
      return level2_rooted_at_r(root, V, terms_dm);
    }


/*
    std::tuple<double, 
               std::unordered_set<int>,
               std::unordered_map<int,int>> 
        level3_alg(const double alpha) {
      // dijktra from the root
      auto p_r = dijkstra(adj, w, root);
      auto dists_r = p_r.first;
      auto trace_r = p_r.second;
      // sort vertices by their distances to the root
      std::vector<int> Lv;
      std::vector<double> Ld;
      for (const auto& p : dists_r) {
        auto v = p.first;
        auto d_rv = p.second;
        if (terms_dm.find(v) != terms_dm.end())
          continue;
        Lv.push_back(v);
        Ld.push_back(d_rv);
      }
      auto idxs_r = argsort(Ld);
      std::vector<int> Lv_sorted;
      std::vector<double> Ld_sorted;
      std::unordered_map<int,size_t> v2i;
      for (size_t i=0; i<idxs_r.size(); i++) {
        auto idx = idxs_r[i];
        Lv_sorted.push_back(Lv[idx]);
        Ld_sorted.push_back(Ld[idx]);
        v2i[Lv[idx]] = i;
      }

      // iteratively add 3-level partial trees
      std::unordered_map<int,int> trace; // trace "short-cut" chosen 2-level trees
      std::unordered_set<int> V_cov, terms_left(terms_dm.begin(), terms_dm.end());
      while (terms_left.size() > 0) {
        // build a 2-level solution as upper bound
        auto best = level2_rooted_at_r(root, V, terms_left);
        if(best.terms_cov.size() == 0)
          break;
        double den2 = best.cost_sc / best.terms_cov.size();

        // enum all u as the single child of the root
        int u_best {NONVERTEX};
        for (auto u: V) {
          // prune by d(r,u)
          if (alpha * best.density() <= den2 - dists_r[u])
            continue;
          
          // collect u's children v by d(u,v)
          std::unordered_set<int> V_; // keep promising v's here
          auto i_u = v2i[u];
          auto d_ru = Ld_sorted[i_u];
          // sweep v's ranked before u
          auto i_v = i_u - 1;
          while (i_v >= 0 and i_u >= 1) {
            auto d_rv = Ld_sorted[i_v];
            if (alpha * d_rv > d_ru - d_rv)
              break;
            V_.insert(Lv_sorted[i_v]);
            i_v -= 1;
          }
          // sweep v's ranked after u
          i_v = i_u + 1;
          while (i_v < Lv.size()) {
            auto d_rv = Ld_sorted[i_v];
            if (alpha * d_rv > d_ru - d_rv) 
              break;
            V_.insert(Lv_sorted[i_v]);
            i_v += 1;
          }

          // add 2-level trees rooted at u 
          double cost_u {d_ru};
          std::unordered_map<int,int> trace_u;
          std::unordered_set<int> terms_left_u(terms_left.begin(), terms_left.end());
          while (terms_left_u.size() > 0) {
            // TODO: will u connect the same v's again?
            auto tree_u = level2_rooted_at_r(u, V_, terms_left_u);
            if (tree_u.terms_cov.size() == 0
                or best.density() <= tree_u.cost_sc/tree_u.terms_cov.size())
              break;
            if (has_key(V_cov, u) and denbest > cost2_u/terms2_u.size()) {
              denbest = cost2_u/terms2_u.size();
              terms_cov_best = std::move(terms2_u);
              trace_best = std::move(trace2_u);
            }

            // TODO: overlap arcs b/w cost_u and cost2_u
            double den_new = (cost_u + cost2_u) / 
              (terms_left.size() - terms_left_u.size() + terms2_u.size());
            if (terms_left.size() > terms_left_u.size()) {
              double den_u = cost_u / (terms_left.size() - terms_left_u.size());
              if (den_u <= den_new)
                break;
            }
            denbest = std::min(denbest, den_new);

            // another LB using terms_left and den(T_u)

            // merge the greedy 2-level partial tree
            cost_u += cost2_u;
            for (auto t: terms2_u) {
              terms_left_u.erase(t);
            }
            for (auto &p: trace2_u) {
              trace_u[p.first] = p.second;
            }
          }

          if (terms_left.size() == terms_left_u.size())
            continue;
          double den_u = cost_u / (terms_left.size() - terms_left_u.size());
          if (den_u < den_min) {
            u_best = u;
            den_min = den_u;
            terms_left_best = std::move(terms_left_u);
            trace_best = std::move(trace_u);
          }
        }

        trace[u_best] = root;
        for (auto &p: trace_best) {
          trace[p.first] = p.second;
        }
        terms_left = std::move(terms_left_best);
      }

      double total {0};
      for (auto e: edges_marked) {
        total += w[{e.first, e.second}];
      }
      std::unordered_set<int> terms_cov;
      std::set_difference(terms_dm.begin(), terms_dm.end(),
                          terms_left.begin(), terms_left.end(),
                          std::inserter(terms_cov, terms_cov.begin()));
      return std::make_tuple(total, terms_cov, trace);
      std::unordered_map<int,int> trace;
      return std::make_tuple(0, std::unordered_set<int> {}, trace);
    }
    */


    double naive_alg() {
      auto p_ = dijkstra(adj, w, root);
      auto dists = p_.first;
      auto trace = p_.second;

      // take the union of all paths from root to t's
      std::set<std::pair<int,int>> edges_marked;
      for (auto t: terms_dm) {
        if (not has_key(trace, t))
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

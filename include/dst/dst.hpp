#pragma once

#include <cmath>
#include <limits>
#include <utility>
#include <cassert>
#include <algorithm>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <tuple>

#include <fmt/ranges.h>
#include <dst/consts.hpp>
#include <dst/utils.hpp>
#include <dst/tree.hpp>


namespace dst {

  std::pair<std::unordered_map<int,double>, 
            std::unordered_map<int,int>> dijkstra(const std::map<int, std::vector<int>> &adj, 
                                                  const std::map<std::pair<int,int>, double> &edgeweight, 
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

      if (has_key(distances, u) and d_u > distances.at(u)) 
        continue;

      trace[u] = u_prev;
      if (not has_key(adj, u))
        continue;
      for (const auto& v: adj.at(u)) {
          double weight = reverse? edgeweight.at({v,u}) : edgeweight.at({u,v});

          if (not has_key(distances, v) or distances.at(v) > distances.at(u) + weight) {
              distances[v] = distances.at(u) + weight;
              pq.emplace(distances.at(v), u, v);
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

      PartialTree::w = &w;
    }


    PartialTree level2_through_v(
        int r,
        int v, 
        double d_rv, 
        const std::unordered_map<int,int> &trace_r,
        const std::unordered_map<int, std::unordered_map<int,double>> &dists_t,
        const std::unordered_map<int, std::unordered_map<int,int>> &trace_t,
        const std::unordered_set<int> &terms_left
    ){
      // collect all t's for the current v
      std::vector<int> terms_left_v;
      std::vector<double> ds_vt;
      for (auto t: terms_left) {
        if (not has_key(dists_t.at(t), v)) 
          continue; // t is not reachable from v
        
        terms_left_v.push_back(t);
        ds_vt.push_back(dists_t.at(t).at(v));
      }

      // sort distances from v to t's
      PartialTree tree {std::make_pair(r, v), d_rv, trace_r};
      auto &&idxs = argsort(ds_vt);
      for (size_t i=0; i < ds_vt.size(); i++) {
        auto idx = idxs[i];
        auto t = terms_left_v[idx];
        double den_i = (tree.cost_sc + ds_vt[idx]) / (i+1);
        if (i == 0 or den_i < tree.density()) {
          PartialTree tree_i {std::make_pair(v, t), 
                              ds_vt[idx],
                              trace_t.at(t),
                              true,
                              true};
          tree.append(tree_i);
        }
        else {
          break; // stop once adding a t increases den_v
        }
      }

      return tree;
    }


    PartialTree level2_rooted_at_r(
        int r,
        const std::unordered_set<int> &V_cand, 
        const std::unordered_set<int> &terms_cand
    ) {
      // dijktra from the root
      auto &&p_r = dijkstra(adj, w, r);
      auto dists_r = std::move(p_r.first);
      auto trace_r = std::move(p_r.second);

      // dijktra from each terminal
      std::unordered_map<int, std::unordered_map<int,int>> trace_t;
      std::unordered_map<int, std::unordered_map<int,double>> dists_t;
      for (auto t: terms_cand) {
        auto &&p_t = dijkstra(adj_r, w, t, true);
        dists_t[t] = std::move(p_t.first);
        trace_t[t] = std::move(p_t.second);
      }

      // iteratively add 2-level partial trees
      PartialTree par {r};
      std::unordered_set<int> terms_left(terms_cand.begin(), terms_cand.end());
      while (terms_left.size() > 0) {
        // enum all v as the middle vertex in a 2-level tree
        PartialTree best {r};
        for (auto v: V_cand) {
          if (not has_key(dists_r, v) or r == v)
            continue;
          //fmt::println("level2_rooted_at_{} v: {}", r, v);
          double d_rv {dists_r.at(v)};
          auto &&tree_v = level2_through_v(r, v, d_rv, trace_r, 
              dists_t, trace_t, terms_left);

          // keep the best across all v
          if (tree_v.density() < best.density())
            best = std::move(tree_v);
        }

        // the rest terminals are not reachable
        if (best.terms_cov.size() == 0)
          break;

        // merge the best 2-level partial tree
        for (auto t: best.terms_cov)
          terms_left.erase(t);
        par.append(best);
      }

      return par;
    }


    PartialTree level2_alg() {
      return level2_rooted_at_r(root, V, terms_dm);
    }


    PartialTree level3_alg(const double alpha) {
      // dijktra from the root
      auto &&p_r = dijkstra(adj, w, root);
      auto dists_r = std::move(p_r.first);
      auto trace_r = std::move(p_r.second);
      // sort vertices by their distances to the root
      std::vector<int> Lv;
      std::vector<double> Ld;
      for (const auto& p : dists_r) {
        auto v = p.first;
        auto d_rv = p.second;
        if (has_key(terms_dm, v))
          continue;
        Lv.push_back(v);
        Ld.push_back(d_rv);
      }
      auto &&idxs_r = argsort(Ld);
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
      std::unordered_set<int> V_cov {root}, 
          terms_left(terms_dm.begin(), terms_dm.end());
      PartialTree tree3 {root};
      while (terms_left.size() > 0) {
        //fmt::println("terms_left: {}", terms_left);
        // build a 2-level solution as upper bound
        auto &&tree2 = level2_rooted_at_r(root, V, terms_left);
        if(tree2.terms_cov.size() == 0)
          break;
        double den2 = tree2.cost_sc / tree2.terms_cov.size();

        // enum all u as the single child of the root
        PartialTree best = std::move(tree2);
        for (auto u: V) {
          // prune by d(r,u)
          if (not has_key(dists_r, u))
            continue;
          // prune by d(r,u)
          if (alpha * best.density() <= den2 - dists_r.at(u))
            continue;
          //fmt::println("u: {}", u);

          // collect u's children v by d(u,v)
          std::unordered_set<int> V_; // keep promising v's here
          auto i_u = v2i.at(u);
          auto d_ru = Ld_sorted[i_u];
          // sweep v's ranked before u
          auto i_v = i_u - 1;
          while (i_v >= 0 and i_u >= 1) {
            auto d_rv = Ld_sorted[i_v];
            if (alpha * d_rv <= d_ru - d_rv)
              break;
            V_.insert(Lv_sorted[i_v]);
            i_v -= 1;
          }
          // sweep v's ranked after u
          i_v = i_u + 1;
          while (i_v < Lv.size()) {
            auto d_rv = Ld_sorted[i_v];
            if (alpha * d_rv <= d_rv - d_ru) 
              break;
            V_.insert(Lv_sorted[i_v]);
            i_v += 1;
          }

          // add 2-level trees rooted at u 
          //fmt::println("V_: {}", V_);
          PartialTree tree_u {std::make_pair(root, u), d_ru, trace_r};
          std::unordered_set<int> terms_left_u(terms_left.begin(), terms_left.end());
          while (terms_left_u.size() > 0) {
            auto &&tree2_u = level2_rooted_at_r(u, V_, terms_left_u);
            if (tree2_u.terms_cov.size() == 0
                or best.density() <= tree2_u.density())
              break;
            if (has_key(V_cov, u) 
                and best.density() > tree2_u.density()) {
              best = tree2_u; // copy
            }

            double den_new = (tree_u.cost_sc + tree2_u.cost_sc) / 
              (terms_left.size() - terms_left_u.size() + tree2_u.terms_cov.size());
            if (terms_left.size() > terms_left_u.size()) { // skip 1st iteration
              double den_u = tree_u.cost_sc / (terms_left.size() - terms_left_u.size());
              if (den_u <= den_new)
                break;
            }

            // TODO: another LB using terms_left and den(T_u)

            // merge tree2_u
            tree_u.append(tree2_u);
            assert(std::abs(den_new - tree_u.density()) < 1e-9);
            if (best.density() > tree_u.density()) {
              best = tree_u; // copy
            }
            for (auto t: tree2_u.terms_cov) {
              terms_left_u.erase(t);
            }
          }
        }

        //fmt::println("best terms_cov: {}", best.terms_cov);
        for (auto t: best.terms_cov)
          terms_left.erase(t);
        tree3.append(best);
      }

      return tree3;
    }


    double naive_alg() {
      auto &&p_ = dijkstra(adj, w, root);
      auto dists = std::move(p_.first);
      auto trace = std::move(p_.second);

      // take the union of all paths from root to t's
      std::set<std::pair<int,int>> edges_marked;
      for (auto t: terms_dm) {
        if (not has_key(trace, t))
            continue; // disconnected graph

        int v {t};
        while (trace.at(v) != NONVERTEX) {
          edges_marked.insert({trace.at(v), v});
          v = trace.at(v);
        }
      }

      double total {0};
      for (auto e: edges_marked)
        total = total + w.at({e.first, e.second});
      return total;
    }

  };

}  // namespace dst

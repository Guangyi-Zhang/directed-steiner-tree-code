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

      if (has_key(distances, u) and d_u > distances.at(u)) 
        continue;

      trace[u] = u_prev;
      for (const auto& v: adj[u]) {
          double weight = reverse? edgeweight.at({v,u}) : edgeweight.at({u,v});

          if (not has_key(distances, v) or distances.at(v) > distances.at(u) + weight) {
              distances[v] = distances.at(u) + weight;
              pq.emplace(distances.at(v), u, v);
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
    int root;
    double cost, cost_sc;
    std::unordered_map<int,int> trace_sc; // trace "short-cut" chosen 2-level trees
    std::unordered_map<int,int> trace;
    std::unordered_set<int> terms_cov;

    PartialTree(int root) : 
        root {root} {
      trace[root] = NONVERTEX;
      trace_sc[root] = NONVERTEX;
    }

    PartialTree(int root, int u, double d_ru) : 
        root {root} {
      trace[root] = NONVERTEX;
      trace_sc[root] = NONVERTEX;
      trace[u] = root;
      trace_sc[u] = root;
      cost_sc += d_ru;
      cost += d_ru;
    }

    void append(PartialTree tree, 
                std::map<std::pair<int,int>, double> edgeweights) {
      // either append r-u-{v}-{t}, or u-{v}-{r}
      cost_sc += tree.cost_sc; // ok to count (u,v) multi times
      auto es = edges();
      for (auto &e: tree.edges()) {
        if (not has_key(es, e)) {
          cost += edgeweights.at(e);
        }
      }

      for (auto t: tree.terms_cov) {
        terms_cov.insert(t);
      }

      for (auto &p: tree.trace_sc) {
        // won't happen: u->v and w->v
        // may happen: -1->u in tree while r->u
        if (not has_key(trace_sc, p.first))
          trace_sc[p.first] = p.second;
      }
      for (auto &p: tree.trace) {
        if (not has_key(trace, p.first))
          trace[p.first] = p.second;
      }
    }

    double density() const {
      if (terms_cov.size() == 0)
        return std::numeric_limits<int>::max();
      return cost_sc / terms_cov.size();
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
        if (not has_key(dists_t.at(t), v)) 
          continue; // t is not reachable from v
        
        terms_left_v.push_back(t);
        ds_vt.push_back(dists_t.at(t).at(v));
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
      PartialTree par {r};
      std::unordered_set<int> terms_left(terms_cand.begin(), terms_cand.end());
      while (terms_left.size() > 0) {
        // used to keep track of the best v so far
        int v_min {NONVERTEX};
        std::vector<int> terms_min;
        double den_min = std::numeric_limits<double>::max();

        // enum all v as the middle vertex in a 2-level tree
        for (auto v: V_cand) {
          if (not has_key(dists_r, v))
            continue;
          double d_rv {dists_r.at(v)};
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
        par.trace_sc[v_min] = r;
        par.cost_sc += dists_r.at(v_min); // ok to count (r,v_min) multi times
        for (auto t: terms_min) {
          par.trace_sc[t] = v_min;
          par.cost_sc += dists_t.at(t).at(v_min);
          terms_left.erase(t);
        }
        // mark the path from root to v
        auto u = v_min;
        while (trace_r.at(u) != NONVERTEX) {
          par.trace[u] = trace_r.at(u);
          u = trace_r.at(u);
        }
        // mark the path from v to each t
        for (auto t: terms_min) {
          u = v_min;
          while (trace_t.at(t).at(u) != NONVERTEX) {
            par.trace[trace_t.at(t).at(u)] = u;
            u = trace_t.at(t).at(u);
          }
        }
      }

      for (auto t: terms_cand) {
        if (not has_key(terms_left, t))
          par.terms_cov.insert(t);
      }
      double total {0};
      for (auto e: par.edges()) {
        total += w.at({e.first, e.second});
      }
      par.cost = total;
      return par;
    }


    PartialTree level2_alg() {
      return level2_rooted_at_r(root, V, terms_dm);
    }


    PartialTree level3_alg(const double alpha) {
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
        if (has_key(terms_dm, v))
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
      std::unordered_set<int> V_cov {root}, 
          terms_left(terms_dm.begin(), terms_dm.end());
      PartialTree tree3 {root};
      while (terms_left.size() > 0) {
        // build a 2-level solution as upper bound
        auto tree2 = level2_rooted_at_r(root, V, terms_left);
        if(tree2.terms_cov.size() == 0)
          break;
        double den2 = tree2.cost_sc / tree2.terms_cov.size();

        // enum all u as the single child of the root
        PartialTree best {tree2}; // copy
        for (auto u: V) {
          // prune by d(r,u)
          if (not has_key(dists_r, u))
            continue;
          if (alpha * best.density() <= den2 - dists_r.at(u))
            continue;
          
          // collect u's children v by d(u,v)
          std::unordered_set<int> V_; // keep promising v's here
          auto i_u = v2i.at(u);
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
          PartialTree tree_u(root, u, d_ru);
          std::unordered_set<int> terms_left_u(terms_left.begin(), terms_left.end());
          while (terms_left_u.size() > 0) {
            auto tree2_u = level2_rooted_at_r(u, V_, terms_left_u);
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
            tree_u.append(tree2_u, w);
            assert(den_new == tree_u.cost_sc);
            if (best.density() > tree_u.density()) {
              best = tree_u; // copy
            }
            for (auto t: tree2_u.terms_cov) {
              terms_left_u.erase(t);
            }
          }
        }

        tree3.append(best, w);
      }

      return tree3;
    }


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

}  // namespace greeter

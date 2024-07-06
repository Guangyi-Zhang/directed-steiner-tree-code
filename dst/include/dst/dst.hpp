#pragma once

#include <cmath>
#include <limits>
#include <memory>
#include <utility>
#include <cassert>
#include <algorithm>
#include <iterator> // std::advance
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <tuple>

#include <fmt/ranges.h>
#include <boost/container_hash/hash.hpp>
#include <dst/consts.hpp>
#include <dst/utils.hpp>
#include <dst/tree.hpp>
#include <dst/tree2.hpp>
#include <dst/dijkstra.hpp>


namespace dst {

  class DST {
  public:
    std::string name;

    int root;
    std::vector<int> terms;
    std::unordered_set<int> terms_dm;// dummy
    std::unordered_map<int, int> terms_map;
    std::unordered_map<std::pair<int,int>, double, boost::hash<std::pair<int,int>>> w;
    std::unordered_map<int, std::vector<int>> adj;
    std::unordered_map<int, std::vector<int>> adj_r; // reverse adj
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
        const std::shared_ptr<std::unordered_map<int,int>> trace_r,
        const std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,double>>>> dists_t,
        const std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>> trace_t,
        const std::unordered_set<int> &terms_left
    ){
      // collect all t's for the current v
      std::vector<int> terms_left_v;
      std::vector<double> ds_vt;
      for (auto t: terms_left) {
        if (not has_key(*(dists_t->at(t)), v)) 
          continue; // t is not reachable from v
        
        terms_left_v.push_back(t);
        ds_vt.push_back(dists_t->at(t)->at(v));
      }

      // sort distances from v to t's
      PartialTree tree {r};
      tree.add_arc(std::make_pair(r, v), d_rv, trace_r);
      auto &&idxs = argsort(ds_vt);
      for (size_t i=0; i < ds_vt.size(); i++) {
        auto idx = idxs[i];
        auto t = terms_left_v[idx];
        double den_i = (tree.cost_sc + ds_vt[idx]) / (i+1);
        if (i == 0 or leq(den_i, tree.density())) {
          tree.add_arc(std::make_pair(v, t), ds_vt[idx], trace_t->at(t), true, true);
        }
        else 
          break; // stop once adding a t increases den_v
      }

      return tree;
    }


    PartialTree level2_rooted_at_r_co(
        int r,
        const std::unordered_set<int> &V_cand, 
        const std::unordered_set<int> &terms_cand
    ) {
      // dijktra from the root
      auto [dists_r, trace_r] = dijkstra(adj, w, r);

      // init a Level2PartialTree for each v
      std::unordered_map<int, Level2PartialTree> trees;
      auto cmp_level2 = [](Level2PartialTree* a, Level2PartialTree* b) { 
        return a->density() < b->density(); 
      };
      std::set<Level2PartialTree*, decltype(cmp_level2)> LBs(cmp_level2); // a set as balanced binary-tree to sort LBs
      for (auto v: V_cand) {
        if (not has_key(*dists_r, v))
          continue;
        trees[v] = std::move(Level2PartialTree {root, v, dists_r->at(v)});
        LBs.insert(&trees.at(v));
      }

      // define what to do after finding a greedy partial tree
      PartialTree par {r};
      int v_best {NONVERTEX};
      CoordinatedDijkstra cosssp {adj_r, w, terms_cand, true}; // dijktra from each terminal
      auto add_greedy = [&] (Level2PartialTree& tree_best) {
        if (DEBUG) fmt::println("level2 greedy through {}, d_rv={}, cov={}, density={}", tree_best.v, tree_best.d_rv, tree_best.terms, tree_best.density());
        PartialTree &&best = tree_best.to_tree(trace_r, cosssp.trace_t);
        for (auto t: best.terms_cov) {
          cosssp.delete_source(t);
        }
        for (auto v: V_cand) {
          if (not has_key(*dists_r, v))
            continue;
          auto &tree_v = trees.at(v);
          tree_v.erase_and_reset(best.terms_cov);
          LBs.insert(&tree_v);
        }
        par.append(best);
        v_best = NONVERTEX; // reset
      };

      // iteratively add 2-level greedy partial trees
      while (terms_cand.size() > par.terms_cov.size()) {
        auto [t, v, d_vt] = cosssp.next();
        if (t == NONVERTEX) { // run out of next()
          for (auto &p: trees) {
            auto &tr = p.second;
            if (v_best == NONVERTEX or 
                lq(tr.density(), trees.at(v_best).density()) or
                (eq(tr.density(), trees.at(v_best).density()) and 
                 tr.terms.size() > trees.at(v_best).terms.size()))
              v_best = tr.v;
          }
          auto &tree_best = trees.at(v_best);
          if (tree_best.terms.size() == 0)
            break;
          add_greedy(tree_best);
        }
        if (has_key(terms_cand, v) or not has_key(*dists_r, v))
          continue;
        auto &tree_v = trees.at(v);
        tree_v.add_term(t, d_vt);

        // update v_best as an UB
        if (v_best == NONVERTEX)
          v_best = v;
        auto &tree_best_old = trees.at(v_best);
        if(tree_v.density() < tree_best_old.density() or
           (eq(tree_v.density(), tree_best_old.density()) and // break ties
            tree_v.terms.size() > tree_best_old.terms.size())) {
          v_best = v;

          Level2PartialTree tmp {tree_v.density()}; // a fake tree
          for (auto it = LBs.upper_bound(&tmp); it != LBs.end(); ) {
            it = LBs.erase(it);
          }
        }
        auto &tree_best = trees.at(v_best);

        // update LB and compare with UB
        if (leq(tree_best.density(), tree_v.density_LB(terms_cand.size() - par.terms_cov.size())))
          LBs.erase(&tree_v);
        else {
          LBs.erase(&tree_v);
          LBs.insert(&tree_v);
        }
        if (LBs.size() >= 2)
          continue;
        if (LBs.size() == 1 and (*LBs.begin())->v != v_best)
          continue;

        // found a greedy partial tree
        if (not tree_best.is_ready(d_vt) and  // until full construction
            tree_best.terms.size() < (terms_cand.size() - par.terms_cov.size()))
          // 1. wait for full construction
          // 2. v_best has reached all remaining terminals
          // 3. some terminals are not reachable
          // 4. v_best need not wait with large d_vt as a LB
          // 5. TODO: bi-direction SSSP
          continue;
        if (tree_best.terms.size() == 0) // the rest terminals not reachable
          break;
        add_greedy(tree_best);
      }

      int sssp_nodes_visited = dists_r->size();
      for (auto &p: *(cosssp.distances_t)) {
        sssp_nodes_visited += p.second->size();
      }
      par.debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      return par;
    }


    PartialTree level2_co_alg() {
      return level2_rooted_at_r_co(root, V, terms_dm);
    }


    PartialTree level2_rooted_at_r(
        int r,
        const std::unordered_set<int> &V_cand, 
        const std::unordered_set<int> &terms_cand,
        std::shared_ptr<std::unordered_map<int,double>> dists_r=nullptr,
        std::shared_ptr<std::unordered_map<int,int>> trace_r=nullptr,
        std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,double>>>> dists_t=nullptr,
        std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>> trace_t=nullptr
    ) {
      // dijktra from the root
      if (dists_r == nullptr)
        std::tie(dists_r, trace_r) = dijkstra(adj, w, r);

      // dijktra from each terminal
      if (dists_t == nullptr) {
        trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
        dists_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,double>>>>();
        for (auto t: terms_cand) {
          auto [dists_, trace_] = dijkstra(adj_r, w, t, true);
          (*dists_t)[t] = dists_; 
          (*trace_t)[t] = trace_;
        }
      }

      // iteratively add 2-level partial trees
      PartialTree par {r};
      int v_best_ {NONVERTEX};
      std::unordered_set<int> terms_left(terms_cand.begin(), terms_cand.end());
      while (terms_left.size() > 0) {
        // enum all v as the middle vertex in a 2-level tree
        PartialTree best {r};
        for (auto v: V_cand) {
          if (not has_key(*dists_r, v) or r == v)
            continue;
          double d_rv {dists_r->at(v)};
          auto &&tree_v = level2_through_v(r, v, d_rv, trace_r, 
              dists_t, trace_t, terms_left);

          // keep the best across all v
          //if (DEBUG) fmt::println("level2_rooted_at_{}: #terms_left={}, trying v={} and density={} and cov={}", r, terms_left.size(), v, tree_v.density(), tree_v.terms_cov);
          if (tree_v.density() < best.density() or // breaking tie
              (std::abs(tree_v.density() - best.density()) < EPSILON and 
               tree_v.terms_cov.size() > best.terms_cov.size())) { //
            best = std::move(tree_v);
            v_best_ = v;
          }
        }

        // the rest terminals are not reachable
        if (best.terms_cov.size() == 0)
          break;

        // merge the best 2-level partial tree
        if (DEBUG) fmt::println("level2_rooted_at_{}: #terms_left={}, best v={} and density={} and cov={}", r, terms_left.size(), v_best_, best.density(), best.terms_cov);
        for (auto t: best.terms_cov)
          terms_left.erase(t);
        par.append(best);
      }

      int sssp_nodes_visited = dists_r->size();
      for (auto &p: *dists_t) {
        sssp_nodes_visited += p.second->size();
      }
      par.debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      return par;
    }


    PartialTree level2_alg() {
      return level2_rooted_at_r(root, V, terms_dm);
    }


    PartialTree level3_alg_outdated(const double alpha) {
      // dijktra from the root
      auto [dists_r, trace_r] = dijkstra(adj, w, root);
      // sort vertices by their distances to the root
      std::vector<int> Lv;
      std::vector<double> Ld;
      for (const auto& p : *dists_r) {
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
      std::unordered_map<int,int> v2i;
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
          if (not has_key(*dists_r, u))
            continue;
          // prune by d(r,u)
          if (alpha * best.density() <= den2 - dists_r->at(u))
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
          PartialTree tree_u {root};
          tree_u.add_arc(std::make_pair(root, u), d_ru, trace_r);
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


    PartialTree naive_alg() {
      auto [dists, trace] = dijkstra(adj, w, root, false, terms_dm);

      // take the union of all paths from root to t's
      PartialTree tree {root};
      for (auto t: terms_dm) {
        if (not has_key(*trace, t))
            continue; // disconnected graph

        tree.add_arc(std::make_pair(root, t), (*dists)[t], trace, false, true);
      }

      tree.debuginfo["sssp_nodes_visited"] = std::to_string(dists->size());
      return tree;
    }

  };

}  // namespace dst

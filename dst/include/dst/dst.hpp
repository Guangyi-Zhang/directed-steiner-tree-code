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
#include <dst/partree.hpp>
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

      Tree::w = &w;
    }


    auto level2_through_v(
        int r,
        int v, 
        double d_rv, 
        const std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,double>>>> dists_t,
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
      auto tree = std::make_shared<PartialTree> (r, v, d_rv);
      auto &&idxs = argsort(ds_vt);
      for (size_t i=0; i < ds_vt.size(); i++) {
        auto idx = idxs[i];
        auto t = terms_left_v[idx];
        double den_i = (tree->cost_sc + ds_vt[idx]) / (i+1);
        if (i == 0 or leq(den_i, tree->density())) {
          tree->add_term(t, ds_vt[idx]);
        }
        else 
          break; // stop once adding a t increases den_v
      }

      return tree;
    }


    auto level2_rooted_at_r_co(
        int r,
        const std::unordered_set<int> &V_cand, 
        const std::unordered_set<int> &terms_cand
    ) {
      // dijktra from the root
      auto [dists_r, trace_r] = dijkstra(adj, w, r);

      // init a PartialTree for each v
      std::unordered_map<int, std::shared_ptr<PartialTree>> trees;
      auto cmp_level2 = [](std::shared_ptr<PartialTree> a, std::shared_ptr<PartialTree> b) { 
        return a->density() < b->density(); 
      };
      std::set<std::shared_ptr<PartialTree>, decltype(cmp_level2)> LBs(cmp_level2); // a set as balanced binary-tree to sort LBs
      for (auto v: V_cand) {
        if (not has_key(*dists_r, v))
          continue;
        trees[v] = std::make_shared<PartialTree> (r, v, dists_r->at(v));
        LBs.insert(trees.at(v));
      }

      // define what to do after finding a greedy partial tree
      auto par = std::make_shared<PartialTreeManager> (r);
      int v_best {NONVERTEX};
      CoordinatedDijkstra cosssp {adj_r, w, terms_cand, true}; // dijktra from each terminal
      auto add_greedy = [&] (std::shared_ptr<PartialTree> best_) {
        auto best = std::make_shared<PartialTree>();
        *best = *best_; // copy
        if (DEBUG) fmt::println("level2 greedy through {}, d_rv={}, cov={}, density={}", best->v, best->d_rv, best->terms, best->density());
        for (auto t: best->terms) {
          cosssp.delete_source(t);
        }
        for (auto v: V_cand) {
          if (not has_key(*dists_r, v))
            continue;
          auto tree_v = trees.at(v);
          tree_v->erase_and_reset(best->terms);
          LBs.insert(tree_v);
        }
        par->append(best);
        v_best = NONVERTEX; // reset
      };

      // iteratively add 2-level greedy partial trees
      while (terms_cand.size() > par->terms.size()) {
        auto [t, v, d_vt] = cosssp.next();
        if (t == NONVERTEX) { // run out of next()
          for (auto &p: trees) {
            auto tr = p.second;
            if (v_best == NONVERTEX or 
                lq(tr->density(), trees.at(v_best)->density()) or
                (eq(tr->density(), trees.at(v_best)->density()) and 
                 tr->terms.size() > trees.at(v_best)->terms.size()))
              v_best = tr->v;
          }
          auto tree_best = trees.at(v_best);
          if (tree_best->terms.size() == 0)
            break;
          add_greedy(tree_best);
        }
        if (has_key(terms_cand, v) or not has_key(*dists_r, v))
          continue;
        auto tree_v = trees.at(v);
        tree_v->add_term(t, d_vt);

        // update v_best as an UB
        if (v_best == NONVERTEX)
          v_best = v;
        auto tree_best_old = trees.at(v_best);
        if(tree_v->density() < tree_best_old->density() or
           (eq(tree_v->density(), tree_best_old->density()) and // break ties
            tree_v->terms.size() > tree_best_old->terms.size())) {
          v_best = v;

          auto tmp = std::make_shared<PartialTree> (tree_v->density()); // a fake tree
          for (auto it = LBs.upper_bound(tmp); it != LBs.end(); ) {
            it = LBs.erase(it);
          }
        }
        auto tree_best = trees.at(v_best);

        // update LB and compare with UB
        if (leq(tree_best->density(), tree_v->density_LB(terms_cand.size() - par->terms.size())))
          LBs.erase(tree_v);
        else {
          LBs.erase(tree_v);
          LBs.insert(tree_v);
        }
        if (LBs.size() >= 2)
          continue;
        if (LBs.size() == 1 and (*LBs.begin())->v != v_best)
          continue;

        // found a greedy partial tree
        if (not tree_best->is_ready(d_vt) and  // until full construction
            tree_best->terms.size() < (terms_cand.size() - par->terms.size()))
          // 1. wait for full construction
          // 2. v_best has reached all remaining terminals
          // 3. some terminals are not reachable
          // 4. v_best need not wait with large d_vt as a LB
          // 5. TODO: bi-direction SSSP
          continue;
        if (tree_best->terms.size() == 0) // the rest terminals not reachable
          break;
        add_greedy(tree_best);
      }

      int sssp_nodes_visited = dists_r->size();
      for (auto &p: *(cosssp.distances_t)) {
        sssp_nodes_visited += p.second->size();
      }
      par->debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      par->trace_r = trace_r;
      par->trace_t = cosssp.trace_t;
      return par;
    }


    auto level2_co_alg() {
      return level2_rooted_at_r_co(root, V, terms_dm);
    }


    auto level2_rooted_at_r(
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
      auto par = std::make_shared<PartialTreeManager> (r);
      int v_best_ {NONVERTEX};
      std::unordered_set<int> terms_left(terms_cand.begin(), terms_cand.end());
      while (terms_left.size() > 0) {
        // enum all v as the middle vertex in a 2-level tree
        std::shared_ptr<PartialTree> best = nullptr;
        for (auto v: V_cand) {
          if (not has_key(*dists_r, v) or r == v)
            continue;
          double d_rv {dists_r->at(v)};
          auto tree_v = level2_through_v(r, v, d_rv, dists_t, terms_left);

          // keep the best across all v
          //if (DEBUG) fmt::println("level2_rooted_at_{}: #terms_left={}, trying v={} and density={} and cov={}", r, terms_left.size(), v, tree_v.density(), tree_v.terms_cov);
          if (best == nullptr or
              tree_v->density() < best->density() or // breaking tie
              (std::abs(tree_v->density() - best->density()) < EPSILON and 
               tree_v->terms.size() > best->terms.size())) { //
            best = tree_v;
            v_best_ = v;
          }
        }

        // the rest terminals are not reachable
        if (best == nullptr or best->terms.size() == 0)
          break;

        // merge the best 2-level partial tree
        if (DEBUG) fmt::println("level2_rooted_at_{}: #terms_left={}, best v={} and density={} and cov={}", r, terms_left.size(), v_best_, best->density(), best->terms);
        for (auto t: best->terms)
          terms_left.erase(t);
        par->append(best);
      }

      int sssp_nodes_visited = dists_r->size();
      for (auto &p: *dists_t) {
        sssp_nodes_visited += p.second->size();
      }
      par->debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      par->trace_r = trace_r;
      par->trace_t = trace_t;
      return par;
    }


    auto level2_alg() {
      return level2_rooted_at_r(root, V, terms_dm);
    }


    auto level3_alg() {
      // pre-compute dijkstra
      // dijktra from the root
      auto [dists_r, trace_r] = dijkstra(adj, w, root);
      // dijkstra from each u
      std::unordered_map<int, Dijkstra> sssp_u;
      for (auto u: V) {
        if (not has_key(*dists_r, u))
          continue;
        sssp_u[u] = std::move(Dijkstra(adj, w, u));
      }
      // backward dijktra from each terminal
      auto trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      auto dists_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,double>>>>();
      for (auto t: terms_dm) {
        if (not has_key(*dists_r, t))
          continue;
        auto [dists_, trace_] = dijkstra(adj_r, w, t, true);
        (*dists_t)[t] = dists_; 
        (*trace_t)[t] = trace_;
      }

      // pre-compute low-level look-up table, one for each v in r-u-v-{t}
      std::unordered_map<int, PartialTreeTable> tbls;
      for (auto v: V) {
        if (v == root)
          continue;
        PartialTreeTable tbl (v);
        for (auto t: terms_dm) {
          auto dists = dists_t->at(t);
          if (not has_key(*dists, v))
            continue;
          tbl.add_term(t, dists->at(v));
        }
        tbl.build();
        tbls[v] = std::move(tbl);
      }

      std::unordered_set<int> terms_left(terms_dm.begin(), terms_dm.end());
      auto tree3 = std::make_shared<PartialTreeManager> (root);
      // iterative add 3-level partial trees
      while (terms_left.size() > 0) {
        std::shared_ptr<PartialTreeManager> best = nullptr;

        for (auto u: V) {
          if (u == root or not has_key(*dists_r, u))
            continue;
          auto &sssp = sssp_u.at(u);
          auto tree_u = std::make_shared<PartialTreeManager> (root, u, dists_r->at(u));
          std::unordered_set<int> terms_left_u(terms_left.begin(), terms_left.end());

          // iterative add 2-level partial trees
          while (terms_left_u.size() > 0) {
            std::shared_ptr<PartialTree> tree2_u = nullptr;
            for (auto v: V) {
              if (v == u or v == root)
                continue;
              while (not has_key(*(sssp.distances), v)) {
                // even after reaching all v, most pq has huge amount left
                auto [z, d_uz] = sssp.next();
                if (z == NONVERTEX)
                  break;
              }
              if (not has_key(*(sssp.distances), v))
                continue;

              double den = tbls.at(v).density(sssp.distances->at(v)); // a lower bound of true tree2_uv
              if (tree2_u == nullptr or tree2_u->density() > den) {
                auto tree2_uv = level2_through_v(u, v, sssp.distances->at(v), dists_t, terms_left_u);
                tree2_u = tree2_uv;
              }
            }
            if (tree2_u == nullptr or tree2_u->terms.size() == 0)
              break;

            double den_new = (tree_u->cost_sc + tree2_u->cost_sc) / 
              (terms_left.size() - terms_left_u.size() + tree2_u->terms.size());
            if (terms_left.size() > terms_left_u.size()) { // skip 1st iteration
              double den_u = tree_u->cost_sc / (terms_left.size() - terms_left_u.size());
              if (den_u <= den_new)
                break;
            }

            // merge tree2_u
            tree_u->append(tree2_u);
            if (best == nullptr or best->density() > tree_u->density()) {
              best = tree_u;
            }
            for (auto t: tree2_u->terms) {
              terms_left_u.erase(t);
            }
          }
        }

        if (best == nullptr or best->terms.size() == 0)
          break;
        for (auto t: best->terms)
          terms_left.erase(t);
        for (auto &p: tbls) {
          p.second.erase(best->terms);
        }
        tree3->append(best);
      }

      int sssp_nodes_visited = dists_r->size();
      for (auto &p: *dists_t) {
        sssp_nodes_visited += p.second->size();
      }
      auto trace_u = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      for (auto &p: sssp_u) {
        sssp_nodes_visited += p.second.distances->size();
        (*trace_u)[p.first] = p.second.trace;
      }
      tree3->debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      tree3->trace_r = trace_r;
      tree3->trace_u = trace_u;
      tree3->trace_t = trace_t;
      return tree3;
    }


    auto level3_alg_naive() {
      // pre-compute dijkstra
      // dijktra from the root
      auto [dists_r, trace_r] = dijkstra(adj, w, root);
      // dijkstra from each u
      auto trace_u = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      auto dists_u = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,double>>>>();
      for (auto u: V) {
        if (u == root or not has_key(*dists_r, u))
          continue;
        auto [dists_, trace_] = dijkstra(adj, w, u);
        (*dists_u)[u] = dists_; 
        (*trace_u)[u] = trace_;
      }
      // backward dijktra from each terminal
      auto trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      auto dists_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,double>>>>();
      for (auto t: terms_dm) {
        if (not has_key(*dists_r, t))
          continue;
        auto [dists_, trace_] = dijkstra(adj_r, w, t, true);
        (*dists_t)[t] = dists_; 
        (*trace_t)[t] = trace_;
      }

      std::unordered_set<int> terms_left(terms_dm.begin(), terms_dm.end());
      auto tree3 = std::make_shared<PartialTreeManager> (root);
      // iterative add 3-level partial trees
      while (terms_left.size() > 0) {
        std::shared_ptr<PartialTreeManager> best = nullptr;

        for (auto u: V) {
          if (u == root or not has_key(*dists_r, u))
            continue;
          auto dists_uv = dists_u->at(u);
          auto tree_u = std::make_shared<PartialTreeManager> (root, u, dists_r->at(u));
          std::unordered_set<int> terms_left_u(terms_left.begin(), terms_left.end());

          // iterative add 2-level partial trees
          while (terms_left_u.size() > 0) {
            std::shared_ptr<PartialTree> tree2_u = nullptr;
            for (auto v: V) {
              if (v == u or v == root or not has_key(*dists_uv, v))
                continue;
              auto tree2_uv = level2_through_v(u, v, dists_uv->at(v), dists_t, terms_left_u);
              if (tree2_u == nullptr or tree2_u->density() > tree2_uv->density())
                tree2_u = tree2_uv;
            }
            if (tree2_u == nullptr or tree2_u->terms.size() == 0)
              break;

            double den_new = (tree_u->cost_sc + tree2_u->cost_sc) / 
              (terms_left.size() - terms_left_u.size() + tree2_u->terms.size());
            if (terms_left.size() > terms_left_u.size()) { // skip 1st iteration
              double den_u = tree_u->cost_sc / (terms_left.size() - terms_left_u.size());
              if (den_u <= den_new)
                break;
            }

            // merge tree2_u
            tree_u->append(tree2_u);
            if (best == nullptr or best->density() > tree_u->density()) {
              best = tree_u;
            }
            for (auto t: tree2_u->terms) {
              terms_left_u.erase(t);
            }
          }
        }

        if (best == nullptr or best->terms.size() == 0)
          break;
        for (auto t: best->terms)
          terms_left.erase(t);
        tree3->append(best);
      }

      int sssp_nodes_visited = dists_r->size();
      for (auto &p: *dists_t) {
        sssp_nodes_visited += p.second->size();
      }
      for (auto &p: *dists_u) {
        sssp_nodes_visited += p.second->size();
      }
      tree3->debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      tree3->trace_r = trace_r;
      tree3->trace_u = trace_u;
      tree3->trace_t = trace_t;
      return tree3;
    }


    auto naive_alg() {
      auto [dists, trace] = dijkstra(adj, w, root, false, terms_dm);

      // take the union of all paths from root to t's
      auto tree = std::make_shared<Tree> (root);
      for (auto t: terms_dm) {
        if (not has_key(*trace, t))
            continue; // disconnected graph

        tree->add_arc(std::make_pair(root, t), (*dists)[t], trace, false, true);
      }

      tree->debuginfo["sssp_nodes_visited"] = std::to_string(dists->size());
      return tree;
    }

  };

}  // namespace dst

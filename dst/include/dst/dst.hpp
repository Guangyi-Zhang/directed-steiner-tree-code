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
#include <iostream>

#include <fmt/ranges.h>
#include <boost/container_hash/hash.hpp>
#include "dst/consts.hpp"
#include "dst/utils.hpp"
#include "dst/tree.hpp"
#include "dst/partree.hpp"
#include "dst/dijkstra.hpp"


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

    // for testing
    int narcs; 
    std::string strtest = "";

    DST(  const std::vector<std::pair<int,int>> &edges, 
          const std::vector<double> &edgeweights, 
          const int root, 
          const std::vector<int> &terminals) :
        root {root}
    {
      int v_max = init_graph(edges, edgeweights);

      // collect unreachble vertices
      auto [dists_r, trace_r] = dijkstra(adj, w, root);

      std::unordered_set<int> unreachables;
      for (auto u: V) {
        if (not has_key(*dists_r, u))
          unreachables.insert(u);
      }

      // re-initialize after proprecessing
      v_max = init_graph(edges, edgeweights, &unreachables);

      for (auto t: terminals) {
        if (has_key(unreachables, t) or not has_key(V, t))
          continue;
        terms.push_back(t);
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


    int init_graph(const std::vector<std::pair<int,int>> &edges, 
                   const std::vector<double> &edgeweights,
                   const std::unordered_set<int> *unreachables=nullptr) 
    {
      V.clear();
      adj.clear();
      adj_r.clear();
      w.clear();
      narcs = 0;

      int v_max = 0;
      for (size_t i = 0; i < edges.size(); i++) {
          auto p = edges[i];

          // drop unreachable arcs
          if (unreachables != nullptr) {
            if (has_key(*unreachables, p.first) or has_key(*unreachables, p.second))
              continue;
          }

          narcs++;
          adj[p.first].push_back(p.second);
          adj_r[p.second].push_back(p.first);
          w[{p.first, p.second}] = edgeweights[i];
          V.insert(p.first);
          V.insert(p.second);
          v_max = std::max({p.first, p.second, v_max});
      }

      return v_max;
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


    auto fast_level2_alg() {
      // dijktra from the root
      auto [dists_r, trace_r] = dijkstra(adj, w, root);

      // init a PartialTree for each v
      std::unordered_map<int, std::shared_ptr<PartialTree>> trees;
      std::list<std::shared_ptr<PartialTree>> LBs;
      for (auto v: V) {
        trees[v] = std::make_shared<PartialTree> (root, v, dists_r->at(v));
        LBs.push_back(trees.at(v));
      }

      CoordinatedDijkstra cosssp {adj_r, w, terms_dm, true}; // dijktra from each terminal
      auto par = std::make_shared<PartialTreeManager> (root);
      par->trace_r = trace_r;
      par->trace_t = cosssp.trace_t;
      int v_best {NONVERTEX};
      auto covered = std::make_shared<std::unordered_set<int>>();

      /*----------------------------------------------
       | define what to do after finding a greedy partial tree
       ----------------------------------------------*/
      auto add_greedy = [&] (std::shared_ptr<PartialTree> best_) {
        auto best = std::make_shared<PartialTree>();
        *best = *best_; // copy
        if (DEBUG) fmt::println("level2 greedy through {}, d_rv={}, cov={}, density={}", best->v, best->d_rv, best->terms, best->density());
        for (auto t: best->terms) {
          cosssp.delete_source(t);
        }
        auto covered_by_best = par->append(best, true);

        // update trees and LBs
        LBs.clear();
        double denbest_next = std::numeric_limits<double>::max();
        for (auto v_tree: trees) { 
          auto tree_v = v_tree.second;
          tree_v->erase_and_reset(best->terms);
          if (has_key(*covered_by_best, tree_v->v) and not has_key(*covered, tree_v->v)) // newly covered
            tree_v->zero_drv();

          if (denbest_next > tree_v->density()) {
            denbest_next = tree_v->density();
            v_best = tree_v->v;
          }

          LBs.push_back(tree_v);
        }

        for (auto v: *covered_by_best)
          covered->insert(v);
      };

      /*----------------------------------------------
       | iteratively add 2-level greedy partial trees
       ----------------------------------------------*/
      while (terms_dm.size() > par->terms.size()) {
        auto [t, v, d_vt] = cosssp.next();
        //fmt::println("t={}, v={}, d_vt={}", t, v, d_vt);

        // run out of next()
        if (t == NONVERTEX) { 
          for (auto &p: trees) { // find best among trees
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
          continue;
        }

        if (has_key(terms_dm, v))
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
        }
        auto tree_best = trees.at(v_best);
        if (DEBUG) fmt::println("best partree so far v={}, cost_sc={}, terms={}, LB={}", v_best, tree_best->cost_sc, tree_best->terms, LBs.front()->density_LB(terms_dm.size() - par->terms.size(), d_vt));

        // compare UB with LBs, remove as many LBs as possible
        while (not LBs.empty() and
               (LBs.front()->ready or 
               leq(tree_best->density(), LBs.front()->density_LB(terms_dm.size() - par->terms.size(), d_vt)))
        ) {
          if (DEBUG) fmt::println("prune v={}, density={}, terms={}", LBs.front()->v, LBs.front()->density(), LBs.front()->terms);
          if (DEBUG_fast_level2) strtest += std::to_string(LBs.front()->v) + ":" + std::to_string(LBs.front()->density_LB(terms_dm.size() - par->terms.size(), d_vt)) + ",";
          LBs.pop_front();
        }
        if (LBs.size() >= 2 or
            (LBs.size() == 1 and LBs.front()->v != v_best))
          continue;

        // found a greedy partial tree
        if (not tree_best->is_ready(d_vt) and  // until full construction
            tree_best->terms.size() < (terms_dm.size() - par->terms.size()))
          // 1. wait for full construction
          // 2. v_best has reached all remaining terminals
          // 3. some terminals are not reachable
          // 4. v_best need not wait with large d_vt as a LB
          // 5. TODO: bi-direction SSSP
          continue;
        if (tree_best->terms.size() == 0) // the rest terminals not reachable
          break;
        if (DEBUG_fast_level2) strtest += "greedy_by_pruning";
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
      par->trace_r = trace_r;
      par->trace_t = trace_t;
      auto covered = std::make_shared<std::unordered_set<int>>();
      std::unordered_set<int> terms_left(terms_cand.begin(), terms_cand.end());
      while (terms_left.size() > 0) {
        // enum all v as the middle vertex in a 2-level tree
        std::shared_ptr<PartialTree> best = nullptr;
        for (auto v: V_cand) {
          double d_rv = has_key(*covered, v)? 0 : dists_r->at(v);
          auto tree_v = level2_through_v(r, v, d_rv, dists_t, terms_left);

          // keep the best across all v
          if (DEBUG) fmt::println("level2_rooted_at_{}: #terms_left={}, trying v={} and density={} and cov={}", r, terms_left.size(), v, tree_v->density(), tree_v->terms);
          if (best == nullptr or
              tree_v->density() < best->density() or // breaking tie
              (std::abs(tree_v->density() - best->density()) < EPSILON and 
               tree_v->terms.size() > best->terms.size())) {
            best = tree_v;
          }
        }

        // the rest terminals are not reachable
        if (best == nullptr or best->terms.size() == 0)
          break;

        // merge the best 2-level partial tree
        if (DEBUG) fmt::println("level2_rooted_at_{}: #terms_left={}, best v={} and density={} and cov={}", r, terms_left.size(), best->v, best->density(), best->terms);
        for (auto t: best->terms)
          terms_left.erase(t);
        auto covered_by_best = par->append(best, true);
        for (auto v: *covered_by_best)
          covered->insert(v);
      }

      int sssp_nodes_visited = dists_r->size();
      for (auto &p: *dists_t) {
        sssp_nodes_visited += p.second->size();
      }
      par->debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      return par;
    }


    auto level2_alg() {
      return level2_rooted_at_r(root, V, terms_dm);
    }


    auto fast_level3_alg(double alpha=1.0, int n_thresholds=10) {
      /*---------------
       | Pre-compute
       --------------*/
      // pre-compute dijkstra
      // dijktra from the root
      auto [dists_r, trace_r] = dijkstra(adj, w, root);
      // dijkstra from each u
      std::unordered_map<int, Dijkstra> sssp_u;
      for (auto u: V) {
        sssp_u[u] = std::move(Dijkstra(adj, w, u));
      }
      // backward dijktra from each terminal
      auto trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      auto dists_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,double>>>>();
      for (auto t: terms_dm) {
        auto [dists_, trace_] = dijkstra(adj_r, w, t, true);
        (*dists_t)[t] = dists_; 
        (*trace_t)[t] = trace_;
      }

      // A pq for each u to rank 2-level partrees
      std::unordered_map<int, 
                         std::priority_queue<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>, 
                                             std::vector<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>>, 
                                             std::greater<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>>>
                        > Q;
      for(auto u: V) {
        Q[u] = std::priority_queue<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>, 
                                     std::vector<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>>, 
                                     std::greater<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>>> ();
      }

      // pre-compute low-level look-up table, one for each v in r-u-v-{t}
      std::unordered_map<int, PartialTreeTable> tbls;
      for (auto v: V) {
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

      // pre-compute min-den among 2-level partrees for each d_uv threshold
      double dist_max = 0; // dist_max may not >= max_{u,v} d_uv
      for (auto &p: *dists_r) {
        dist_max = std::max(dist_max, p.second);
      }
      auto thrminden = ThresholdedMinDensity(n_thresholds, dist_max, tbls);

      /*---------------------------------
       | Start constructing 3-level tree
       --------------------------------*/
      std::unordered_set<int> terms_left(terms_dm.begin(), terms_dm.end());
      auto tree3 = std::make_shared<PartialTreeManager> (root);
      auto trace_u = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      for (auto &p: sssp_u) {
        (*trace_u)[p.first] = p.second.trace;
      }
      tree3->trace_r = trace_r;
      tree3->trace_u = trace_u;
      tree3->trace_t = trace_t;
      auto covered = std::make_shared<std::unordered_set<int>>();
      // ***** iterative add 3-level partial trees *****
      while (terms_left.size() > 0) {
        // best 2-level partree as initial UB 
        auto best = std::make_shared<PartialTreeManager> (root);
        std::shared_ptr<PartialTree> best2 = nullptr;
        for (auto v: V) {
          double d_rv = has_key(*covered, v)? 0 : dists_r->at(v);
          auto tree2_rv = level2_through_v(root, v, d_rv, dists_t, terms_left);
          if (best2 == nullptr or lq(tree2_rv->density(), best2->density()))
            best2 = tree2_rv;
        }
        if (best2 == nullptr or best2->terms.size() == 0)
          break;
        if (DEBUG) fmt::println("best2 v={}, cost_sc={}, #terms={}", best2->v, best2->cost_sc, best2->terms.size());
        best->append(best2);

        // ***** find the best 3-level partree *****
        for (auto u: V) {
          double d_ru = has_key(*covered, u)? 0 : dists_r->at(u);
          if (leq(alpha * best->density(), best2->density() - d_ru)) {
            // it is fine to skip root=u, same as best2
            // also fine to skip u in covered with d_ru=0 
            continue;
          }

          auto Q_u = &(Q.at(u));
          auto &sssp = sssp_u.at(u);
          auto tree_u = std::make_shared<PartialTreeManager> (root, u, d_ru);
          std::shared_ptr<PartialTree> tree2_u_old = nullptr;
          std::unordered_set<int> terms_left_u;
          for (auto t: terms_left) {
            if (has_key(*(dists_t->at(t)), u))
              terms_left_u.insert(t);
          }
          int niter_u = 0;
          bool is_u_pruned = false;

          // ***** iterative add 2-level partial trees *****
          while (terms_left_u.size() > 0) {
            niter_u++;
            std::shared_ptr<PartialTree> tree2_u = nullptr;

            // ***** process prior reached v's *****
            if (not Q_u->empty() and tree2_u_old == nullptr) { // restore Q_u at the i-th visit to u 
              auto Q_new = std::priority_queue<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>, 
                                               std::vector<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>>, 
                                               std::greater<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>>> ();

              std::vector<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>> &wrap_Q_u = Container(*Q_u);
              for(auto it = wrap_Q_u.begin(); it != wrap_Q_u.end(); it++) { // iterate Q_u
                auto [den, niter, d_uv, v, tree2_uv] = *it;
                tree2_uv = tbls.at(v).partree(u, d_uv, &terms_left_u);
                Q_new.emplace(tree2_uv->density(), niter_u, d_uv, v, tree2_uv);
              }
              Q[u] = std::move(Q_new); Q_u = &(Q.at(u));

              auto [den, niter, d_uv, v, tree2_uv] = Q_u->top();
              if (tree2_u == nullptr or lq(tree2_uv->density(), tree2_u->density())) {
                tree2_u = tree2_uv;
              }
            }
            else if (not Q_u->empty()){ // update Q_u across iterations within the i-th visit
              auto [den, niter, d_uv, v, tree2_uv] = Q_u->top();
              while(niter != niter_u) {
                Q_u->pop();
                if (tree2_uv == nullptr) {
                  tree2_uv = tbls.at(v).partree(u, d_uv, &terms_left_u);
                } else {
                  tree2_uv->erase_and_reset(tree2_u_old->terms);
                }
                Q_u->emplace(tree2_uv->density(), niter_u, d_uv, v, tree2_uv);
                std::tie(den, niter, d_uv, v, tree2_uv) = Q_u->top();
              }

              if (tree2_u == nullptr or lq(tree2_uv->density(), tree2_u->density())) {
                tree2_u = tree2_uv;
              }
            }

            // ***** process new v's *****
            while (true) { 
              // TODO: setup to_reach in sssp
              auto [v, d_uv] = sssp.next(); // even after reaching all v, most pq has huge amount left
              if (v == NONVERTEX)
                break;
              if (u != root and v == root) // allow u == v
                continue;

              // try to early-terminate sssp
              double denlb = thrminden.min_density(d_uv);
              if (tree2_u != nullptr and leq(tree2_u->density(), denlb))
                break;
              if (leq(alpha * best->density(), denlb)) {
                is_u_pruned = true;
                break;
              }
              double tree_u_LB = (tree_u->cost_sc + denlb * terms_left_u.size()) / (terms_left_u.size() + tree_u->terms.size());
              if (leq(alpha * best->density(), tree_u_LB)) {
                is_u_pruned = true;
                break;
              }

              double den = tbls.at(v).density(d_uv); // a lower bound of true tree2_uv
              if (tree2_u != nullptr and leq(tree2_u->density(), den)) {
                Q_u->emplace(den, niter_u, d_uv, v, nullptr); // lazy evaluate tree2_uv
                continue;
              }

              auto tree2_uv = tbls.at(v).partree(u, d_uv, &terms_left_u);
              Q_u->emplace(tree2_uv->density(), niter_u, d_uv, v, tree2_uv);
              if (tree2_u == nullptr or lq(tree2_uv->density(), tree2_u->density())) {
                tree2_u = tree2_uv;
              }
            }

            // criterion for refusing new tree2_u
            if (is_u_pruned or tree2_u == nullptr or tree2_u->terms.size() == 0)
              break;
            double den_new = (tree_u->cost_sc + tree2_u->cost_sc) / 
              (tree_u->terms.size() + tree2_u->terms.size());
            if (not tree_u->terms.empty() and lq(tree_u->density(), den_new))
              break;

            // ***** merge tree2_u: best among tree2_uv *****
            if (DEBUG) fmt::println("u={} adds partree rooted at v={}, cost_sc={}, #terms={}", u, tree2_u->v, tree2_u->cost_sc, tree2_u->terms.size());
            auto tree2_u_cp = tree2_u->copy();
            tree_u->append(tree2_u_cp);
            if (DEBUG) fmt::println("tree_u cost_sc={}, #terms={}", tree_u->cost_sc, tree_u->terms.size());
            if (best == nullptr or lq(tree_u->density(), best->density())) {
              best = tree_u;
            }
            for (auto t: tree2_u->terms) {
              terms_left_u.erase(t);
            }
            tree2_u_old = tree2_u_cp;
          } // end of processing vertex u
        } // found a greedy 3-level partree

        if (best == nullptr or best->terms.size() == 0)
          break;
        for (auto t: best->terms)
          terms_left.erase(t);
        for (auto &p: tbls) {
          p.second.erase(best->terms);
        }
        thrminden = ThresholdedMinDensity(n_thresholds, dist_max, tbls);
        auto covered_by_best = tree3->append(best, true);
        tree3->density1by1.push_back(best->density()); // for teting
        tree3->costsc1by1.push_back(best->cost_sc); // for testing
        for (auto v: *covered_by_best)
          covered->insert(v);
        if (DEBUG) fmt::println("add partree r-u={}, r-v={}, cost_sc={}, #terms={}", best->u, (best->subtrees.size() > 0)? best->subtrees[0].first->v: -1, best->cost_sc, best->terms.size());
      }

      int sssp_nodes_visited = dists_r->size();
      for (auto &p: *dists_t) {
        sssp_nodes_visited += p.second->size();
      }
      for (auto &p: sssp_u) {
        sssp_nodes_visited += p.second.distances->size();
      }
      tree3->debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      return tree3;
    }


    auto level3_alg() {
      // pre-compute dijkstra
      // dijktra from the root
      auto [dists_r, trace_r] = dijkstra(adj, w, root);
      // dijkstra from each u
      auto trace_u = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      auto dists_u = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,double>>>>();
      for (auto u: V) {
        auto [dists_, trace_] = dijkstra(adj, w, u);
        (*dists_u)[u] = dists_; 
        (*trace_u)[u] = trace_;
      }
      // backward dijktra from each terminal
      auto trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      auto dists_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,double>>>>();
      for (auto t: terms_dm) {
        auto [dists_, trace_] = dijkstra(adj_r, w, t, true);
        (*dists_t)[t] = dists_; 
        (*trace_t)[t] = trace_;
      }

      std::unordered_set<int> terms_left(terms_dm.begin(), terms_dm.end());
      auto tree3 = std::make_shared<PartialTreeManager> (root);
      tree3->trace_r = trace_r;
      tree3->trace_u = trace_u;
      tree3->trace_t = trace_t;
      auto covered = std::make_shared<std::unordered_set<int>>();
      // iterative add 3-level partial trees
      while (terms_left.size() > 0) {
        std::shared_ptr<PartialTreeManager> best = nullptr;

        // loop each u
        for (auto u: V) {
          auto dists_uv = dists_u->at(u);
          double d_ru = has_key(*covered, u)? 0 : dists_r->at(u);
          auto tree_u = std::make_shared<PartialTreeManager> (root, u, d_ru);
          std::unordered_set<int> terms_left_u;
          for (auto t: terms_left) {
            if (has_key(*(dists_t->at(t)), u))
              terms_left_u.insert(t);
          }

          // iterative add 2-level partial trees by looping each v
          while (terms_left_u.size() > 0) {
            std::shared_ptr<PartialTree> tree2_u = nullptr;
            for (auto v: V) {
              if ((u != root and v == root) or not has_key(*dists_uv, v)) // allow u == v
                continue;
              auto tree2_uv = level2_through_v(u, v, dists_uv->at(v), dists_t, terms_left_u);
              if (tree2_u == nullptr or tree2_u->density() > tree2_uv->density())
                tree2_u = tree2_uv;
            } // end of looping v
            if (tree2_u == nullptr or tree2_u->terms.size() == 0)
              break;

            // try the best tree2_u
            double den_new = (tree_u->cost_sc + tree2_u->cost_sc) / 
              (tree_u->terms.size() + tree2_u->terms.size());
            if (not tree_u->terms.empty() and lq(tree_u->density(), den_new))
              break;

            // merge tree2_u
            tree_u->append(tree2_u);
            if (best == nullptr or lq(tree_u->density(), best->density())) {
              best = tree_u;
            }
            for (auto t: tree2_u->terms) {
              terms_left_u.erase(t);
            }
          } // end of one u
        } // end of looping u

        if (best == nullptr or best->terms.size() == 0)
          break;
        for (auto t: best->terms)
          terms_left.erase(t);
        auto covered_by_best = tree3->append(best, true);
        tree3->density1by1.push_back(best->density()); // for teting
        tree3->costsc1by1.push_back(best->cost_sc); // for testing
        for (auto v: *covered_by_best)
          covered->insert(v);
        if (DEBUG) fmt::println("add partree r-u={}, cost_sc={}, #terms={}", best->u, best->cost_sc, best->terms.size());
      }

      int sssp_nodes_visited = dists_r->size();
      for (auto &p: *dists_t) {
        sssp_nodes_visited += p.second->size();
      }
      for (auto &p: *dists_u) {
        sssp_nodes_visited += p.second->size();
      }
      tree3->debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      return tree3;
    }


    auto adaptive_level1_alg() {

      auto tree = std::make_shared<Tree> (root);
      int sssp_nodes_visited = 0;

      std::unordered_set<int> terms_left(terms_dm.begin(), terms_dm.end());
      while (not terms_left.empty()) {
        auto [dists, trace] = dijkstra(adj, w, root, false, terms_left);
        sssp_nodes_visited += dists->size();

        // find the best shortest path
        int t_min = NONVERTEX;
        double dist_min = std::numeric_limits<double>::max();
        for (auto t: terms_left) {
          if (not has_key(*dists, t))
            continue; // disconnected graph
          if (dist_min > dists->at(t)) {
            dist_min = dists->at(t);
            t_min = t;
          }
        }

        // include the best shortest path
        if (t_min != NONVERTEX) {
          tree->add_arc(std::make_pair(root, t_min), (*dists)[t_min], trace, nullptr, false, true);
          terms_left.erase(t_min);

          // zero weight for covered arcs
          auto v = t_min;
          while (trace->at(v) != NONVERTEX) {
            auto e = std::make_pair(trace->at(v), v);
            w[e] = 0;
            v = trace->at(v);
          }
        } 
        else {
          break;
        }
      }

      tree->debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      return tree;
    }


    auto level1_alg() {
      auto [dists, trace] = dijkstra(adj, w, root, false, terms_dm);

      // take the union of all paths from root to t's
      auto tree = std::make_shared<Tree> (root);
      for (auto t: terms_dm) {
        if (not has_key(*trace, t))
          continue; // disconnected graph

        tree->add_arc(std::make_pair(root, t), (*dists)[t], trace, nullptr, false, true);
      }

      tree->debuginfo["sssp_nodes_visited"] = std::to_string(dists->size());
      return tree;
    }

  };

}  // namespace dst

#pragma once

#include <cmath>
#include <limits>
#include <memory>
#include <utility>
#include <cassert>
#include <algorithm>
#include <iterator> // std::advance
#include <numeric>  // std::partial_sum
#include <string>
#include <vector>
#include <queue>
#include <functional>
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
#include "dst/csrgraph.hpp"


namespace dst {

  class DST {
  public:
    std::string name;

    int root;
    std::vector<int> terms;
    std::unordered_set<int> terms_dm;// dummy
    std::unordered_map<int, int> terms_map;
    CSRGraph graph;
    int n_original = 0;

    // for testing
    int narcs;
    std::string strtest = "";

    DST(  const std::vector<std::pair<int,int>> &edges,
          const std::vector<float> &edgeweights,
          const int root,
          const std::vector<int> &terminals) :
        root {root}
    {
      // First pass: build temporary CSR to find unreachable vertices
      int v_max = build_csr(graph, edges, edgeweights, nullptr, {});

      // Dijkstra from root on the temporary graph
      auto [dists_r, trace_r] = dijkstra(graph.fwd_offset, graph.fwd_target,
                                          graph.fwd_weight, root);

      std::unordered_set<int> unreachables;
      for (int u = 0; u <= v_max; u++) {
        if (not has_key(*dists_r, u))
          unreachables.insert(u);
      }

      // Filter terminals
      for (auto t: terminals) {
        if (has_key(unreachables, t))
          continue;
        if (t < 0 or t > v_max)
          continue;
        terms.push_back(t);
      }

      // Build dummy terminal edges
      std::vector<std::pair<std::pair<int,int>, float>> dummy_edges;
      int i = 1;
      for (auto t: terms) {
        int t_dm = v_max + i++;
        terms_dm.insert(t_dm);
        terms_map[t] = t_dm;
        terms_map[t_dm] = t;
        dummy_edges.push_back({{t, t_dm}, 0.0f});
      }

      // Second pass: build final CSR with filtered edges + dummy terminals
      v_max = build_csr(graph, edges, edgeweights, &unreachables, dummy_edges);
      n_original = v_max + 1 - static_cast<int>(terms.size());
      narcs = graph.m - static_cast<int>(dummy_edges.size());

      Tree::csr_graph = &graph;
    }


    static int build_csr(
        CSRGraph &g,
        const std::vector<std::pair<int,int>> &edges,
        const std::vector<float> &edgeweights,
        const std::unordered_set<int> *unreachables,
        const std::vector<std::pair<std::pair<int,int>, float>> &extra_edges)
    {
      // Pass 1: determine n and count degrees
      int v_max = 0;
      int edge_count = 0;
      for (size_t i = 0; i < edges.size(); i++) {
        auto &p = edges[i];
        if (unreachables != nullptr) {
          if (has_key(*unreachables, p.first) or has_key(*unreachables, p.second))
            continue;
        }
        v_max = std::max({p.first, p.second, v_max});
        edge_count++;
      }
      for (auto &[ep, ew] : extra_edges) {
        v_max = std::max({ep.first, ep.second, v_max});
        edge_count++;
      }

      int n = v_max + 1;
      int m = edge_count;
      g.n = n;
      g.m = m;

      // Allocate
      g.fwd_offset.assign(n + 1, 0);
      g.rev_offset.assign(n + 1, 0);
      g.fwd_target.resize(m);
      g.fwd_weight.resize(m);
      g.rev_target.resize(m);
      g.rev_weight.resize(m);

      // Count out-degree (fwd) and in-degree (rev)
      for (size_t i = 0; i < edges.size(); i++) {
        auto &p = edges[i];
        if (unreachables != nullptr) {
          if (has_key(*unreachables, p.first) or has_key(*unreachables, p.second))
            continue;
        }
        g.fwd_offset[p.first + 1]++;
        g.rev_offset[p.second + 1]++;
      }
      for (auto &[ep, ew] : extra_edges) {
        g.fwd_offset[ep.first + 1]++;
        g.rev_offset[ep.second + 1]++;
      }

      // Prefix sum
      for (int i = 1; i <= n; i++) {
        g.fwd_offset[i] += g.fwd_offset[i-1];
        g.rev_offset[i] += g.rev_offset[i-1];
      }

      // Fill using temporary position counters, just a copy of fwd_offset
      std::vector<int> fwd_pos(g.fwd_offset.begin(), g.fwd_offset.end() - 1);
      std::vector<int> rev_pos(g.rev_offset.begin(), g.rev_offset.end() - 1);

      for (size_t i = 0; i < edges.size(); i++) {
        auto &p = edges[i];
        if (unreachables != nullptr) {
          if (has_key(*unreachables, p.first) or has_key(*unreachables, p.second))
            continue;
        }
        float w = edgeweights[i];
        int fi = fwd_pos[p.first]++;
        g.fwd_target[fi] = p.second;
        g.fwd_weight[fi] = w;
        int ri = rev_pos[p.second]++;
        g.rev_target[ri] = p.first;
        g.rev_weight[ri] = w;
      }
      for (auto &[ep, ew] : extra_edges) {
        int fi = fwd_pos[ep.first]++;
        g.fwd_target[fi] = ep.second;
        g.fwd_weight[fi] = ew;
        int ri = rev_pos[ep.second]++;
        g.rev_target[ri] = ep.first;
        g.rev_weight[ri] = ew;
      }

      return v_max;
    }


    /*
     * level2_through_v builds one 2-level bunch root->v->{t}: terminals reachable from v are
     * added in increasing v-t distance until density would rise.
     */
    auto level2_through_v(
        int r,
        int v,
        double d_rv,
        const std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,float>>>> dists_t,
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
      auto [dists_r, trace_r] = dijkstra(graph.fwd_offset, graph.fwd_target,
                                          graph.fwd_weight, root);

      // init a Level2CenterIndex for each v
      std::unordered_map<int, Level2CenterIndex> center_idxs;
      std::list<int> LBs;
      auto covered_terms = std::make_shared<std::unordered_set<int>>();
      for (int v = 0; v < n_original; v++) {
        if (not has_key(*dists_r, v))
          continue;
        center_idxs.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(v),
          std::forward_as_tuple(root, v, dists_r->at(v), covered_terms)
        );
        LBs.push_back(v);
      }

      CoordinatedDijkstra cosssp {graph.rev_offset, graph.rev_target,
                                   graph.rev_weight, terms_dm};
      auto par = std::make_shared<PartialTreeManager> (root);
      par->trace_r = trace_r;
      par->trace_t = cosssp.trace_t;
      int v_best {NONVERTEX};
      auto covered = std::make_shared<std::unordered_set<int>>();

      /*----------------------------------------------
       | define what to do after finding a greedy partial tree
       ----------------------------------------------*/
      auto add_greedy = [&] (int greedy_v) {
        auto best = std::make_shared<PartialTree>();
        *best = center_idxs.at(greedy_v).tree; // copy
        if (DEBUG) fmt::println("level2 greedy through {}, d_rv={}, cov={}, density={}", best->v, best->d_rv, best->terms, best->density());
        for (auto t: best->terms) {
          cosssp.delete_source(t);
          covered_terms->insert(t);
        }
        auto covered_by_best = par->append(best, true);

        // update center_idxs and LBs
        LBs.clear();
        double denbest_next = std::numeric_limits<double>::max();
        for (auto &[v, idx]: center_idxs) {
          idx.erase_and_reset(best->terms);
          if (has_key(*covered_by_best, idx.tree.v) and not has_key(*covered, idx.tree.v)) // newly covered
            idx.zero_drv();

          if (denbest_next > idx.density()) {
            denbest_next = idx.density();
            v_best = idx.tree.v;
          }

          LBs.push_back(v);
        }

        for (auto v: *covered_by_best)
          covered->insert(v);
      };

      /*----------------------------------------------
       | iteratively add 2-level greedy partial trees
       ----------------------------------------------*/
      while (terms_dm.size() > par->terms.size()) {
        auto [t, v, d_vt] = cosssp.next();

        // run out of next()
        if (t == NONVERTEX) {
          for (auto &[v_, idx]: center_idxs) { // find best among trees
            if (v_best == NONVERTEX or
                lq(idx.density(), center_idxs.at(v_best).density()) or
                (eq(idx.density(), center_idxs.at(v_best).density()) and
                 idx.tree.terms.size() > center_idxs.at(v_best).tree.terms.size()))
              v_best = idx.tree.v;
          }
          if (center_idxs.at(v_best).tree.terms.size() == 0)
            break;
          add_greedy(v_best);
          continue;
        }

        if (has_key(terms_dm, v))
          continue;
        auto &idx_v = center_idxs.at(v);
        idx_v.add_term(t, d_vt);

        // update v_best as an UB
        if (v_best == NONVERTEX)
          v_best = v;
        auto &idx_best_old = center_idxs.at(v_best);
        if(idx_v.density() < idx_best_old.density() or
           (eq(idx_v.density(), idx_best_old.density()) and
            idx_v.tree.terms.size() > idx_best_old.tree.terms.size())) {
          v_best = v;
        }
        auto &idx_best = center_idxs.at(v_best);
        int k = terms_dm.size() - par->terms.size();
        if (DEBUG) fmt::println("best partree so far v={}, cost_sc={}, terms={}, LB={}", v_best, idx_best.tree.cost_sc, idx_best.tree.terms, center_idxs.at(LBs.front()).tree.density_LB(k, d_vt));

        // compare UB with LBs, remove as many LBs as possible
        while (not LBs.empty()) {
          auto &idx_lb = center_idxs.at(LBs.front());
          auto lb = idx_lb.optimistic_bound(idx_lb.tree.d_rv, k, d_vt, false);
          if (idx_lb.tree.ready or leq(idx_best.density(), lb.density)) {
            if (DEBUG) fmt::println("prune v={}, density={}, terms={}", idx_lb.tree.v, idx_lb.density(), idx_lb.tree.terms);
            if (DEBUG_fast_level2) strtest += std::to_string(idx_lb.tree.v) + ":" + std::to_string(lb.density) + ",";
            LBs.pop_front();
          } else {
            break;
          }
        }
        if (LBs.size() >= 2 or
            (LBs.size() == 1 and LBs.front() != v_best))
          continue;

        // found a greedy partial tree
        if (not idx_best.is_ready(d_vt) and  // until full construction
            idx_best.tree.terms.size() < (terms_dm.size() - par->terms.size()))
          continue;
        if (idx_best.tree.terms.size() == 0) // the rest terminals not reachable
          break;
        if (DEBUG_fast_level2) strtest += "greedy_by_pruning";
        #ifndef NO_FAST2_LB
        add_greedy(v_best);
        #endif
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
        int n_V,
        const std::unordered_set<int> &terms_cand,
        std::shared_ptr<std::unordered_map<int,float>> dists_r=nullptr,
        std::shared_ptr<std::unordered_map<int,int>> trace_r=nullptr,
        std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,float>>>> dists_t=nullptr,
        std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>> trace_t=nullptr
    ) {
      // dijktra from the root
      if (dists_r == nullptr)
        std::tie(dists_r, trace_r) = dijkstra(graph.fwd_offset, graph.fwd_target,
                                               graph.fwd_weight, r);

      // dijktra from each terminal
      if (dists_t == nullptr) {
        trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
        dists_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,float>>>>();
        for (auto t: terms_cand) {
          auto [dists_, trace_] = dijkstra(graph.rev_offset, graph.rev_target,
                                            graph.rev_weight, t);
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
        for (int v = 0; v < n_V; v++) {
          if (not has_key(*dists_r, v))
            continue;
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
      return level2_rooted_at_r(root, n_original, terms_dm);
    }


    auto tree_vertices(const std::shared_ptr<Tree> tree) const {
      std::unordered_set<int> vertices;
      for (auto &p: tree->trace) {
        vertices.insert(p.first);
        if (p.second != NONVERTEX)
          vertices.insert(p.second);
      }
      return vertices;
    }


    auto residual_weights(const std::shared_ptr<Tree> tree) const {
      // Set all tree edges to cost 0
      auto fwd_res = graph.fwd_weight;
      auto rev_res = graph.rev_weight;
      for (auto &p: tree->trace) {
        int v = p.first;
        int u = p.second;
        if (u == NONVERTEX)
          continue;
        for (int i = graph.fwd_offset[u]; i < graph.fwd_offset[u+1]; i++) {
          if (graph.fwd_target[i] == v) { fwd_res[i] = 0.0f; break; }
        }
        for (int i = graph.rev_offset[v]; i < graph.rev_offset[v+1]; i++) {
          if (graph.rev_target[i] == u) { rev_res[i] = 0.0f; break; }
        }
      }
      return std::make_pair(std::move(fwd_res), std::move(rev_res));
    }


    auto adaptive_level2_alg() {
      // Terminal-side distances are fixed in the root-adaptive variant.
      auto trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      auto dists_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,float>>>>();
      for (auto t: terms_dm) {
        auto [dists_, trace_] = dijkstra(graph.rev_offset, graph.rev_target,
                                            graph.rev_weight, t);
          (*dists_t)[t] = dists_;
          (*trace_t)[t] = trace_;
        }

        auto par = std::make_shared<PartialTreeManager> (root);
      par->trace_t = trace_t;
      std::unordered_set<int> terms_left(terms_dm.begin(), terms_dm.end());
      int sssp_nodes_visited = 0;
      for (auto &p: *dists_t) {
        sssp_nodes_visited += p.second->size();
      }

      while (terms_left.size() > 0) {
        auto sources = tree_vertices(par->totree);
        auto [dists_r, trace_r] = dijkstra_multisource(graph.fwd_offset, graph.fwd_target,
                                                         graph.fwd_weight, sources);
        sssp_nodes_visited += dists_r->size();
        par->trace_r = trace_r;

        std::shared_ptr<PartialTree> best = nullptr;
        for (int v = 0; v < n_original; v++) {
          if (not has_key(*dists_r, v))
            continue;
          auto tree_v = level2_through_v(root, v, dists_r->at(v), dists_t, terms_left);

          if (best == nullptr or
              tree_v->density() < best->density() or
              (eq(tree_v->density(), best->density()) and
               tree_v->terms.size() > best->terms.size())) {
            best = tree_v;
          }
        }

        if (best == nullptr or best->terms.size() == 0)
          break;

        for (auto t: best->terms)
          terms_left.erase(t);
        par->append(best, true);
      }

      par->debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      return par;
    }


    auto full_adaptive_level2_alg() {
      // Root-adaptive variant (adaptive_level2) fixes terminal-side distances;
      // here both sides are recomputed each iteration on residual weights.
      auto par = std::make_shared<PartialTreeManager> (root);
      std::unordered_set<int> terms_left(terms_dm.begin(), terms_dm.end());
      int sssp_nodes_visited = 0;

      while (terms_left.size() > 0) {
        auto [fwd_res, rev_res] = residual_weights(par->totree); // Paid tree edges cost 0
        auto [dists_r, trace_r] = dijkstra(graph.fwd_offset, graph.fwd_target,
                                          fwd_res, root);

        // Backward Dijkstra from each remaining terminal on residual weights.
        auto trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
        auto dists_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,float>>>>();
        for (auto t: terms_left) {
          auto [dists_, trace_] = dijkstra(graph.rev_offset, graph.rev_target,
                                            rev_res, t);
          (*dists_t)[t] = dists_;
          (*trace_t)[t] = trace_;
        }

        sssp_nodes_visited += dists_r->size();
        for (auto &p: *dists_t) {
          sssp_nodes_visited += p.second->size();
        }
        par->trace_r = trace_r;
        par->trace_t = trace_t;

        // Greedy 2-level bunch: pick v minimizing density among uncovered terms.
        std::shared_ptr<PartialTree> best = nullptr;
        for (int v = 0; v < n_original; v++) {
          if (not has_key(*dists_r, v))
            continue;
          auto tree_v = level2_through_v(root, v, dists_r->at(v), dists_t, terms_left);

          if (best == nullptr or
              tree_v->density() < best->density() or
              (eq(tree_v->density(), best->density()) and
               tree_v->terms.size() > best->terms.size())) {
            best = tree_v;
          }
        }

        if (best == nullptr or best->terms.size() == 0)
          break;

        for (auto t: best->terms)
          terms_left.erase(t);
        par->append(best, true);
      }

      par->debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      return par;
    }


    // One greedy step: pick a 2-level bunch (root->v->terminals) minimizing density.
    auto fast_level2_minden(
        std::shared_ptr<std::unordered_map<int,float>> dists_r,
        const std::unordered_set<int> &terms_left,
        const std::vector<float> &rev_weight_to_use
    ) {
      // init a Level2CenterIndex for each v reachable from root
      std::unordered_map<int, Level2CenterIndex> center_idxs;
      std::list<int> LBs;
      for (int v = 0; v < n_original; v++) {
        if (not has_key(*dists_r, v))
          continue;
        center_idxs.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(v),
          std::forward_as_tuple(root, v, dists_r->at(v))
        );
        LBs.push_back(v);
      }

      CoordinatedDijkstra cosssp {graph.rev_offset, graph.rev_target,
                                   rev_weight_to_use, terms_left};
      int v_best {NONVERTEX};

      // function: update the v_best with lowest density (UB) if tree_v has lower density
      auto update_best = [&] (int v) {
        if (v_best == NONVERTEX) {
          v_best = v;
          return;
        }
        auto &idx_v = center_idxs.at(v);
        auto &idx_best = center_idxs.at(v_best);
        if (idx_v.density() < idx_best.density() or
            (eq(idx_v.density(), idx_best.density()) and
             idx_v.tree.terms.size() > idx_best.tree.terms.size())) {
          v_best = v;
        }
      };

      // function: choose the v_best among all centers with lowest density (UB)
      auto choose_best = [&] () -> std::shared_ptr<PartialTree> {
        for (auto &p: center_idxs) {
          update_best(p.first);
        }
        if (v_best == NONVERTEX)
          return nullptr;
        auto best = std::make_shared<PartialTree>();
        *best = center_idxs.at(v_best).tree; // copy
        return best;
      };

      while (true) {
        auto [t, v, d_vt] = cosssp.next();
        if (t == NONVERTEX) { // coordinated SSSP exhausted
          auto best = choose_best();
          int sssp_nodes_visited = 0;
          for (auto &p: *(cosssp.distances_t))
            sssp_nodes_visited += p.second->size();
          return std::make_tuple(best, cosssp.trace_t, sssp_nodes_visited);
        }

        if (has_key(terms_dm, v) or not has_key(center_idxs, v))
          continue;

        auto &idx_v = center_idxs.at(v);
        idx_v.add_term(t, d_vt);
        update_best(v);
        auto &idx_best = center_idxs.at(v_best);

        // Prune from the queue front only. Later centers do not need to be
        // inspected until every earlier LB has been certified/pruned.
        int k = terms_left.size();
        while (not LBs.empty()) {
          auto &idx_lb = center_idxs.at(LBs.front());
          auto lb = idx_lb.optimistic_bound(idx_lb.tree.d_rv, k, d_vt, false);
          if (idx_lb.tree.ready or leq(idx_best.density(), lb.density)) {
            LBs.pop_front();
          } else {
            break;
          }
        }

        // wait until only v_best remains among surviving LBs
        if (LBs.size() >= 2 or
            (LBs.size() == 1 and LBs.front() != v_best))
          continue;

        // early stop once the best tree is fully determined
        if (not idx_best.is_ready(d_vt) and
            idx_best.tree.terms.size() < terms_left.size())
          continue;

        #ifndef NO_FAST2_LB
        auto best = std::make_shared<PartialTree>();
        *best = idx_best.tree; // copy
        int sssp_nodes_visited = 0;
        for (auto &p: *(cosssp.distances_t))
          sssp_nodes_visited += p.second->size();
        return std::make_tuple(best, cosssp.trace_t, sssp_nodes_visited);
        #endif
      }
    }


    // Root-adaptive fast level-2: fixed terminal-side distances, dynamic root-side
    // distances. One persistent backward CoordinatedDijkstra; per-center Fenwick
    // indexes maintain best terminal prefixes as discoveries arrive.
    auto adaptive_fast_level2_alg() {
      auto par = std::make_shared<PartialTreeManager> (root);
      par->trace_t = nullptr;
      std::unordered_set<int> terms_left(terms_dm.begin(), terms_dm.end());
      if (terms_left.empty()) {
        par->debuginfo["sssp_nodes_visited"] = "0";
        return par;
      }

      // One lazily allocated index per root-reachable center v: terminals are
      // inserted in nondecreasing d(v,t) order.
      std::unordered_map<int, AdaptiveLevel2CenterIndex> center_idxs;

      // Terminal-side search runs once for the whole algorithm run.
      CoordinatedDijkstra cosssp {graph.rev_offset, graph.rev_target,
                                   graph.rev_weight, terms_dm};
      par->trace_t = cosssp.trace_t;
      int root_sssp_nodes_visited = 0;
      double terminal_lb = 0; // frontier distance from cosssp; lower-bounds undiscovered d(v,t)
      bool terminal_search_exhausted = false;

      // Function: prefer lower density; tie-break by covering more terminals.
      auto better_prefix = [](
          double a_density,
          int a_nterm,
          double b_density,
          int b_nterm
      ) {
        if (a_nterm == 0)
          return false;
        if (b_nterm == 0)
          return true;
        return a_density < b_density or
               (eq(a_density, b_density) and a_nterm > b_nterm);
      };

      auto bound_can_beat = [&](
          const Level2PruneBound &bound,
          double best_density,
          int best_nterm
      ) {
        return better_prefix(bound.density, bound.nterm, best_density, best_nterm);
      };

      // Get or create the center index for vertex v, lazily populated from cosssp distances
      auto ensure_center = [&] (int v) -> AdaptiveLevel2CenterIndex& {
        // Try to emplace; `inserted` is true only if v is new to the map
        auto [it, inserted] = center_idxs.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(v),
            std::forward_as_tuple(v, static_cast<int>(terms_dm.size()))
        );
        // Already exists → return existing entry without re-initializing
        if (!inserted)
          return it->second;

        // Collect (distance, terminal) pairs for v from cosssp, filtered to remaining terms
        std::vector<std::pair<double, int>> discovered;
        for (auto &p: *(cosssp.distances_t)) {
          int t = p.first;
          if (not has_key(terms_left, t) or not has_key(*(p.second), v))
            continue;
          discovered.emplace_back(p.second->at(v), t);
        }
        // Sort by distance ascending and populate the new center index
        std::sort(discovered.begin(), discovered.end());
        for (auto &[d_vt, t]: discovered)
          it->second.insert(t, d_vt);

        return it->second;
      };

      while (terms_left.size() > 0) {
        // Root-side residual distances from vertices already in the partial tree.
        auto sources = tree_vertices(par->totree);
        auto [dists_r, trace_r] = dijkstra_multisource(graph.fwd_offset, graph.fwd_target,
                                                         graph.fwd_weight, sources);
        root_sssp_nodes_visited += dists_r->size();
        par->trace_r = trace_r;

        // Survivor list of still-possible centers for this greedy iteration.
        std::list<int> survivors;
        std::unordered_set<int> survivors_set;
        int best_v = NONVERTEX;
        double best_density = std::numeric_limits<double>::max();
        int best_nterm = 0;

        auto refresh_center = [&] (int v) {
          auto cand = center_idxs.at(v).discovered_bound(dists_r->at(v));
          if (better_prefix(cand.density, cand.nterm, best_density, best_nterm)) {
            best_v = v;
            best_density = cand.density;
            best_nterm = cand.nterm;
          }
        };

        // Initialize the survivors list and the best center
        for (auto &p: *dists_r) {
          int v = p.first;
          if (v >= n_original)
            continue;
          auto &idx = ensure_center(v);
          survivors.push_back(v);
          survivors_set.insert(v);
          auto cand = idx.discovered_bound(p.second);
          if (better_prefix(cand.density, cand.nterm, best_density, best_nterm)) {
            best_v = v;
            best_density = cand.density;
            best_nterm = cand.nterm;
          }
        }

        // Advance terminal-side search until the best prefix is certified.
        while (true) {
          if (best_nterm > 0) {
            while (not survivors.empty()) {
              int v = survivors.front();
              auto bound = center_idxs.at(v).optimistic_bound(
                  dists_r->at(v),
                  terms_left.size(),
                  terminal_lb,
                  terminal_search_exhausted
              );
              if (bound_can_beat(bound, best_density, best_nterm))
                break; // not yet certified

              survivors_set.erase(v);
              survivors.pop_front();
            }

            if (survivors.empty()) {
              #ifndef NO_FAST2_LB
              break; // no survivors, certified
              #endif
            }
          }

          if (terminal_search_exhausted)
            break;

          auto [t, v, d_vt] = cosssp.next();
          if (t == NONVERTEX) {
            terminal_search_exhausted = true;
            terminal_lb = std::numeric_limits<double>::infinity();
            continue;
          }

          terminal_lb = d_vt;
          if (v < n_original and has_key(*dists_r, v)) {
            auto &idx = ensure_center(v);
            idx.insert(t, d_vt);
            if (has_key(survivors_set, v)) // v is still a survivor
              refresh_center(v);
          }
        }

        auto best = best_v == NONVERTEX?
            std::shared_ptr<PartialTree> {nullptr} :
            center_idxs.at(best_v).partree(root, dists_r->at(best_v));
        if (best == nullptr or best->terms.size() == 0)
          break;

        // Remove covered terminals from the persistent search and all indexes.
        for (auto t: best->terms) {
          terms_left.erase(t);
          cosssp.delete_source(t);
        }
        for (auto &p: center_idxs)
          p.second.erase(best->terms);
        par->append(best, true);
      }

      int sssp_nodes_visited = root_sssp_nodes_visited;
      for (auto &p: *(cosssp.distances_t)) {
        sssp_nodes_visited += p.second->size();
      }
      par->debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      return par;
    }


    auto full_adaptive_fast_level2_alg() {
      auto par = std::make_shared<PartialTreeManager> (root);
      std::unordered_set<int> terms_left(terms_dm.begin(), terms_dm.end());
      int sssp_nodes_visited = 0;

      while (terms_left.size() > 0) {
        auto [fwd_res, rev_res] = residual_weights(par->totree);
        auto [dists_r, trace_r] = dijkstra(graph.fwd_offset, graph.fwd_target,
                                            fwd_res, root);
        sssp_nodes_visited += dists_r->size();

        auto [best, trace_t, sssp_t] = fast_level2_minden(dists_r, terms_left, rev_res);
        sssp_nodes_visited += sssp_t;
        par->trace_r = trace_r;
        par->trace_t = trace_t;

        if (best == nullptr or best->terms.size() == 0)
          break;

        for (auto t: best->terms)
          terms_left.erase(t);
        par->append(best, true);
      }

      par->debuginfo["sssp_nodes_visited"] = std::to_string(sssp_nodes_visited);
      return par;
    }


    auto fast_level3_alg(double alpha=1.0, int n_thresholds=10) {
      /*---------------
       | Pre-compute
       --------------*/
      // pre-compute dijkstra
      // dijktra from the root
      auto [dists_r, trace_r] = dijkstra(graph.fwd_offset, graph.fwd_target,
                                          graph.fwd_weight, root);
      // dijkstra from each u
      std::unordered_map<int, Dijkstra> sssp_u;
      for (int u = 0; u < n_original; u++) {
        if (not has_key(*dists_r, u))
          continue;
        sssp_u[u] = std::move(Dijkstra(graph.fwd_offset, graph.fwd_target,
                                         graph.fwd_weight, u));
      }
      // backward dijktra from each terminal
      auto trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      auto dists_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,float>>>>();
      for (auto t: terms_dm) {
        auto [dists_, trace_] = dijkstra(graph.rev_offset, graph.rev_target,
                                          graph.rev_weight, t);
        (*dists_t)[t] = dists_;
        (*trace_t)[t] = trace_;
      }

      // A pq for each u to rank 2-level partrees
      std::unordered_map<int,
                         std::priority_queue<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>,
                                             std::vector<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>>,
                                             std::greater<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>>>
                        > Q;
      for (auto &p : sssp_u) {
        Q[p.first] = std::priority_queue<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>,
                                     std::vector<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>>,
                                     std::greater<std::tuple<double,int,double,int,std::shared_ptr<PartialTree>>>> ();
      }

      // pre-compute low-level look-up table, one for each v in r-u-v-{t}
      std::unordered_map<int, PartialTreeTable> tbls;
      for (int v = 0; v < n_original; v++) {
        if (not has_key(*dists_r, v))
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

      // pre-compute min-den among 2-level partrees for each d_uv threshold
      float dist_max = 0; // dist_max may not >= max_{u,v} d_uv
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
        for (int v = 0; v < n_original; v++) {
          if (not has_key(*dists_r, v))
            continue;
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
        for (auto &[u, sssp_entry] : sssp_u) {
          if (not has_key(*dists_r, u))
            continue;
          double d_ru = has_key(*covered, u)? 0 : dists_r->at(u);
          #ifndef NO_FAST3_PRUNE_U
          if (leq(alpha * best->density(), best2->density() - d_ru)) {
            continue;
          }
          #endif

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
              if (has_key(terms_dm, v))
                continue;

              // try to early-terminate sssp
              double denlb = thrminden.min_density(d_uv);
              #ifndef NO_FAST3_LB
              if (tree2_u != nullptr and leq(tree2_u->density(), denlb))
                break;
              #endif
              #ifndef NO_FAST3_PRUNE_SUBTREE
              if (leq(alpha * best->density(), denlb)) {
                is_u_pruned = true;
                break;
              }
              #endif
              #ifndef NO_FAST3_PRUNE_SUBTREE
              double tree_u_LB = (tree_u->cost_sc + denlb * terms_left_u.size()) / (terms_left_u.size() + tree_u->terms.size());
              if (leq(alpha * best->density(), tree_u_LB)) {
                is_u_pruned = true;
                break;
              }
              #endif

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
      auto [dists_r, trace_r] = dijkstra(graph.fwd_offset, graph.fwd_target,
                                          graph.fwd_weight, root);
      // dijkstra from each u
      auto trace_u = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      auto dists_u = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,float>>>>();
      for (int u = 0; u < n_original; u++) {
        if (not has_key(*dists_r, u))
          continue;
        auto [dists_, trace_] = dijkstra(graph.fwd_offset, graph.fwd_target,
                                          graph.fwd_weight, u);
        (*dists_u)[u] = dists_;
        (*trace_u)[u] = trace_;
      }
      // backward dijktra from each terminal
      auto trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      auto dists_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,float>>>>();
      for (auto t: terms_dm) {
        auto [dists_, trace_] = dijkstra(graph.rev_offset, graph.rev_target,
                                          graph.rev_weight, t);
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
        for (int u = 0; u < n_original; u++) {
          if (not has_key(*dists_u, u))
            continue;
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
            for (int v = 0; v < n_original; v++) {
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
        auto [dists, trace] = dijkstra(graph.fwd_offset, graph.fwd_target,
                                        graph.fwd_weight, root, terms_left);
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
            int u = trace->at(v);
            graph.zero_edge(u, v);
            v = u;
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
      auto [dists_r, trace_r] = dijkstra(graph.fwd_offset, graph.fwd_target,
                                          graph.fwd_weight, root, terms_dm);

      // take the union of all paths from root to t's
      auto tree = std::make_shared<Tree> (root);
      for (auto t: terms_dm) {
        if (not has_key(*trace_r, t))
          continue; // disconnected graph

        tree->add_arc(std::make_pair(root, t), (*dists_r)[t], trace_r, nullptr, false, true);
      }

      tree->debuginfo["sssp_nodes_visited"] = std::to_string(dists_r->size());
      return tree;
    }

  };

}  // namespace dst

#pragma once

#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <list>
#include <vector>

#include "dst/utils.hpp"
#include "dst/tree.hpp"


namespace dst {

  enum class Level2BoundMode { Cheap, Threshold };

  struct Level2PrefixBound {
    double density {std::numeric_limits<double>::max()};
    int nterm {0};
  };

  inline double level2_threshold(int q, double next_weight, double prefix_weight) {
    return q * next_weight - prefix_weight;
  }

  template <class ThresholdFn>
  int best_level2_prefix_length(int total_terms, double d_rv, ThresholdFn threshold) {
    if (total_terms <= 1)
      return total_terms;
    int best_q = 0;
    int left = 1, right = total_terms - 1;
    while (left <= right) {
      int mid = left + (right - left) / 2;
      if (leq(threshold(mid), d_rv)) {
        best_q = mid;
        left = mid + 1;
      } else {
        right = mid - 1;
      }
    }
    return best_q + 1;
  }

  class BasePartialTree {
    public:
    int root {NONVERTEX}; 
    std::unordered_set<int> terms;
    double cost_sc = 0;

    std::shared_ptr<std::unordered_map<int,int>> trace_r = nullptr;
    std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>> trace_t = nullptr;
    std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>> trace_u = nullptr;

    std::unordered_map<std::string, std::string> debuginfo;

    BasePartialTree() {};

    double density() {
      if (terms.size() == 0)
        return std::numeric_limits<double>::max();
      return cost_sc / terms.size();
    }
  };


  /*
    2-level partial tree that supports density LB, terminal removal
  */
  class PartialTree: public BasePartialTree {
    public:
    int v {NONVERTEX}; // the root has a single child
    double d_rv;

    std::unordered_map<int, double> distances_t;
    std::list<int> terms_after_ready, terms_before_ready;
    bool ready = false; 

    PartialTree() {};

    PartialTree(int r, int v, double d_rv) : 
        v {v}, d_rv {d_rv} {
      root = r;
      cost_sc = d_rv;
    }

    void add_term(int t, double d_vt) {
      distances_t[t] = d_vt;
      double new_den = (cost_sc + d_vt) / (terms.size() + 1);
      if (ready or lq(density(), new_den)) {
        ready = true;
        terms_after_ready.push_back(t);
      } else {
        terms.insert(t);
        terms_before_ready.push_back(t);
        cost_sc += d_vt;
      }
    }

    void erase_and_reset(const std::unordered_set<int> &terms_del) {
      if (&terms == &terms_del) { // erase itself
        for (auto t: terms_del) {
          cost_sc -= distances_t.at(t);
        }
        terms.clear();
      } else {
        for (auto t: terms_del) {
          if (not has_key(terms, t))
            continue;
          terms.erase(t);
          cost_sc -= distances_t.at(t);
        }
      }

      // reset
      ready = false;
      while (not ready) {
        if (terms_after_ready.empty())
          break;
        auto t = terms_after_ready.front();
        terms_after_ready.pop_front();

        if (has_key(terms_del, t)) 
          continue;

        add_term(t, distances_t.at(t));
        if (ready) {
          terms_after_ready.pop_back();
          terms_after_ready.push_front(t);
        }
      }
    }

    void zero_drv() {
      // reset
      cost_sc = 0;
      ready = false;
      std::unordered_set<int> terms_cp = std::move(terms);
      std::list<int> terms_before_ready_cp = std::move(terms_before_ready);
      std::list<int> terms_after_ready_cp = std::move(terms_after_ready);
      terms.clear();
      terms_before_ready.clear();
      terms_after_ready.clear();

      // re-add items in order
      // no need to iterate t's in terms_after_ready
      std::list<int> terms_excluded;
      for (auto t: terms_before_ready_cp) {
        if (not has_key(terms_cp, t)) 
          continue;

        if (ready)
          terms_excluded.push_front(t);
        else 
          add_term(t, distances_t.at(t));
      }
      terms_after_ready.splice(terms_after_ready.end(), terms_excluded);
      terms_after_ready.splice(terms_after_ready.end(), terms_after_ready_cp);
    }

    double density_LB(int k, double d_vt_lb=0) {
      if (ready)
        return density();
      else
        return (cost_sc + (k-terms.size()) * d_vt_lb) / k;
    }

    bool is_ready(double d_LB=-1) {
      if (ready)
        return true;

      if (d_LB >= 0) {
        double new_den = (cost_sc + d_LB) / (terms.size() + 1);
        if (lq(density(), new_den)) 
          ready = true;
      }

      return ready;
    }

    auto copy() {
      auto cp = std::make_shared<PartialTree>();
      *cp = *this;
      return cp;
    }

    auto to_tree () {
      // TODO: test
      auto tree = std::make_shared<Tree> (root);
      tree->add_arc(std::make_pair(root, v), d_rv, trace_r);
      for (auto t: terms) {
        tree->add_arc(std::make_pair(v, t), distances_t.at(t), trace_t->at(t), nullptr, true, true);
      }
      return tree;
    }
  };


  // Binary indexed tree (BIT): O(log n) point add and prefix sum over a fixed-size array.
  class FenwickTree {
    public:
    std::vector<double> bit;     // internal sums with 1-indexed position
    std::vector<double> values;  // current value at each 0-indexed position

    FenwickTree() {};

    FenwickTree(int n) {
      reset(n);
    }

    void reset(int n) {
      bit.assign(n + 1, 0);
      values.assign(n, 0);
    }

    int size() const {
      return values.size();
    }

    double at(int idx) const {
      assert(0 <= idx and idx < size());
      return values[idx];
    }

    // Add delta to values[idx] and propagate into bit.
    void add(int idx, double delta) {
      assert(0 <= idx and idx < size());
      values[idx] += delta;
      // bit[i] (1-indexed) holds the sum of the last (i & -i) elements ending at i.
      // i & -i is 2^k (lowest set bit), so each hop is a power-of-two jump; binary carry
      // clears low bits and bumps the next, giving at most O(log n) nodes.
      // bit[i] == values[i-1] only when i is odd (block length 1); at i = 2^k it sums 2^k
      // elements ending at i; other i sum intermediate-sized blocks (e.g. i=6 → values[4..5]).
      for (int i = idx + 1; i < static_cast<int>(bit.size()); i += i & -i)
        bit[i] += delta;
    }

    void set(int idx, double value) {
      add(idx, value - at(idx));
    }

    // Sum of values[0..idx].
    double prefix_sum(int idx) const {
      if (idx < 0 or size() == 0)
        return 0;
      if (idx >= size())
        idx = size() - 1;

      double sum = 0;
      // Prefix values[0..idx] = 1-indexed positions [1..idx+1]. Walk right-to-left,
      // peeling off disjoint power-of-two blocks: bit[i] is the sum of the last (i & -i)
      // elements ending at i (the rightmost chunk not yet counted). i -= i & -i skips
      // past that chunk; the chunks partition [1..idx+1] exactly, so we never double-count.
      // e.g. idx=4 → i=5 adds [5], i=4 adds [1..4].
      for (int i = idx + 1; i > 0; i -= i & -i)
        sum += bit[i];
      return sum;
    }

    double total_sum() const {
      return prefix_sum(size() - 1);
    }

    // Smallest index whose prefix sum is >= target (0-indexed rank).
    int select_by_cumulative(double target) const {
      if (target <= 0 or target > total_sum() + EPSILON)
        return NONVERTEX;

      // Largest power-of-two jump that still fits in bit[].
      int step = 1;
      while ((step << 1) < static_cast<int>(bit.size()))
        step <<= 1;

      // Find smallest 0-indexed rank r with prefix_sum(r) >= target (inverse of prefix_sum).
      // Same Fenwick blocks, assembled left-to-right: at each step size, try extending idx
      // by step via bit[idx+step] (sum of the next step elements in block order). Commit
      // the jump if acc + bit[idx+step] < target (k-th item not reached yet); otherwise
      // skip and halve step to search inside the block where target falls.
      // e.g. counts [0,0,1,0,0,1,0,0], target=2: step=4, idx=0 → 0+bit[4]=1 < 2, commit;
      // step=2, idx=4 → 1+bit[6]=2 >= 2, skip and refine with smaller step.
      int idx = 0;
      double acc = 0;
      for (; step > 0; step >>= 1) {
        int next = idx + step;
        if (next < static_cast<int>(bit.size()) and acc + bit[next] < target) {
          idx = next;
          acc += bit[next];
        }
      }

      return idx;
    }
  };


  // Per-center index for fast level-2: wraps PartialTree with closed-form
  // threshold mode. Fenwick trees are not needed here because terminals arrive
  // in nondecreasing d(v,t) order and PartialTree already tracks the best
  // prefix: when ready, the prefix is final; when not ready, all discovered
  // alive terminals are in tree.terms.
  class Level2CenterIndex {
    public:
    PartialTree tree;

    Level2CenterIndex() {};

    Level2CenterIndex(int r, int v, double d_rv) :
        tree {r, v, d_rv} {};

    void add_term(int t, double d_vt) {
      tree.add_term(t, d_vt);
    }

    void erase_and_reset(const std::unordered_set<int> &terms_del) {
      tree.erase_and_reset(terms_del);
    }

    void zero_drv() {
      tree.zero_drv();
    }

    double density() { return tree.density(); }
    bool is_ready(double d_LB = -1) { return tree.is_ready(d_LB); }

    Level2PrefixBound optimistic_bound(
        double /*d_rv*/,
        int total_alive_terms,
        double terminal_lb,
        bool terminal_search_exhausted,
        Level2BoundMode mode
    ) const {
      int k = total_alive_terms;
      if (k == 0)
        return Level2PrefixBound {};

      int m = static_cast<int>(tree.terms.size());
      double C = tree.cost_sc;

      if (tree.ready) {
        Level2PrefixBound bound;
        bound.density = C / m;
        bound.nterm = m;
        return bound;
      }

      // Not ready: all discovered alive terminals are in tree.terms.
      int unknown = terminal_search_exhausted ? 0 : std::max(0, k - m);

      if (m == 0) {
        Level2PrefixBound bound;
        bound.density = (C + k * terminal_lb) / k;
        bound.nterm = k;
        return bound;
      }

      if (unknown == 0) {
        Level2PrefixBound bound;
        bound.density = C / m;
        bound.nterm = m;
        return bound;
      }

      // m > 0, unknown > 0: Cheap fills the tail; Threshold picks the better.
      double full_density = (C + unknown * terminal_lb) / k;

      if (mode == Level2BoundMode::Cheap) {
        Level2PrefixBound bound;
        bound.density = full_density;
        bound.nterm = k;
        return bound;
      }

      // Threshold: full tail beats current prefix iff m * lb <= C.
      double current_density = C / m;
      Level2PrefixBound bound;
      if (leq(full_density, current_density)) {
        bound.density = full_density;
        bound.nterm = k;
      } else {
        bound.density = current_density;
        bound.nterm = m;
      }
      return bound;
    }
  };


  // Per-center index for adaptive level-2: terminals in fixed nondecreasing d(v,t)
  // order, with dynamic alive/covered status via two Fenwick trees.
  class AdaptiveLevel2CenterIndex {
    public:
    int v {NONVERTEX};
    std::vector<int> term_order;              // rank -> terminal id
    std::vector<double> distances;            // rank -> fixed d(v, t)
    std::unordered_map<int, int> rank_of_term; // terminal -> rank
    FenwickTree alive_counts;  // 1 per alive terminal, 0 if covered
    FenwickTree alive_weights; // d(v,t) per alive terminal, 0 if covered

    AdaptiveLevel2CenterIndex() {};

    AdaptiveLevel2CenterIndex(int v, int max_terms) :
        v {v}, alive_counts {max_terms}, alive_weights {max_terms} {};

    // Append a newly discovered terminal at the next rank (nondecreasing d(v,t)).
    void insert(int t, double d_vt) {
      if (has_key(rank_of_term, t))
        return;

      int rank = term_order.size();
      assert(rank < alive_counts.size());
      term_order.push_back(t);
      distances.push_back(d_vt);
      rank_of_term[t] = rank;
      alive_counts.add(rank, 1);
      alive_weights.add(rank, d_vt);
    }

    // Mark covered terminals dead in both Fenwick trees (values stay at rank).
    void erase(const std::unordered_set<int> &terms_del) {
      for (auto t: terms_del) {
        if (not has_key(rank_of_term, t))
          continue;

        int rank = rank_of_term.at(t);
        if (alive_counts.at(rank) <= 0.5) // floating point comparison
          continue;
        alive_counts.add(rank, -1);
        alive_weights.add(rank, -distances[rank]);
      }
    }

    int alive_discovered_count() const {
      return static_cast<int>(alive_counts.total_sum() + 0.5); // floating point rounding
    }

    // Weight sum of the first nterm alive terminals in fixed rank order.
    double alive_prefix_weight(int nterm) const {
      if (nterm <= 0)
        return 0;
      int rank = alive_counts.select_by_cumulative(nterm);
      if (rank == NONVERTEX)
        return 0;
      return alive_weights.prefix_sum(rank);
    }

    double threshold(int q) const {
      int rank_next = alive_counts.select_by_cumulative(q + 1);
      return level2_threshold(q, distances[rank_next], alive_prefix_weight(q));
    }

    int best_prefix_length(double d_rv) const {
      int m = alive_discovered_count();
      auto _threshold = [&](int q) -> double { return threshold(q); };
      return best_level2_prefix_length(m, d_rv, _threshold);
    }

    // Density bound using only discovered alive terminals.
    Level2PrefixBound discovered_bound(double d_rv) const {
      Level2PrefixBound bound;
      int nterm = best_prefix_length(d_rv);
      if (nterm == 0)
        return bound;

      bound.nterm = nterm;
      bound.density = (d_rv + alive_prefix_weight(nterm)) / nterm;
      return bound;
    }

    // Optimistic density lower bound for this center.
    //
    // Discovered alive terminals have exact d(v,t) values in rank order. Any
    // undiscovered alive terminal has not been popped from the global
    // CoordinatedDijkstra frontier yet, so its d(v,t) is at least terminal_lb.
    // For certification we ask: even if every undiscovered alive terminal were
    // exactly terminal_lb, what is the best prefix density this center could
    // still achieve? If that optimistic density cannot beat the current best,
    // this center is safe to discard for the current greedy iteration.
    Level2PrefixBound optimistic_bound(
        double d_rv,
        int total_alive_terms,
        double terminal_lb,
        bool terminal_search_exhausted,
        Level2BoundMode mode = Level2BoundMode::Threshold
    ) const {
      int m = alive_discovered_count();
      int unknown = terminal_search_exhausted ? 0 : std::max(0, total_alive_terms - m);
      int total = m + unknown;
      if (total == 0)
        return Level2PrefixBound {};

      if (mode == Level2BoundMode::Cheap) {
        Level2PrefixBound bound;
        bound.nterm = total;
        bound.density = (d_rv + alive_weights.total_sum() + unknown * terminal_lb) / total;
        return bound;
      }

      auto _threshold = [&](int q) -> double {
        if (q < m) {
          int rank_next = alive_counts.select_by_cumulative(q + 1);
          return level2_threshold(q, distances[rank_next], alive_prefix_weight(q));
        }
        return level2_threshold(m, terminal_lb, alive_weights.total_sum());
      };

      int nterm = best_level2_prefix_length(total, d_rv, _threshold);

      Level2PrefixBound bound;
      bound.nterm = nterm;
      bound.density = (
          d_rv +
          alive_prefix_weight(std::min(nterm, m)) +
          std::max(0, nterm - m) * terminal_lb
      ) / nterm;
      return bound;
    }

    // Build PartialTree(root, v, d_rv) with the best discovered alive prefix.
    auto partree(int root, double d_rv) const {
      auto par = std::make_shared<PartialTree> (root, v, d_rv);
      int nterm = best_prefix_length(d_rv);
      for (int i = 1; i <= nterm; i++) {
        int rank = alive_counts.select_by_cumulative(i);
        int t = term_order[rank];
        double d_vt = distances[rank];
        par->terms.insert(t);
        par->distances_t[t] = d_vt;
        par->cost_sc += d_vt;
      }
      return par;
    }
  };


  class PartialTreeManager: public BasePartialTree {
    // can be used a 3-level partial tree if u != NONVERTEX
    public:
    int u {NONVERTEX}; // the root has a single child
    double d_ru;
    std::shared_ptr<Tree> totree = nullptr;
    // use vector, as trees are processed (add_arc()) by the order they are added
    // bool flag to mark if it has been added into totree
    std::vector<std::pair<std::shared_ptr<PartialTree>, bool>> subtrees; // used as r->{v}->{t}
    std::vector<std::pair<std::shared_ptr<PartialTreeManager>, bool>> trees; // used as r->{u}->{v}->{t}

    std::list<double> density1by1; // for testing
    std::list<double> costsc1by1; // for testing

    PartialTreeManager() {}; // only for nullptr

    PartialTreeManager(int r) {
      root = r;
      totree = std::make_shared<Tree> (root);
    };

    PartialTreeManager(int r, int u, double d_ru) : 
        u {u}, d_ru {d_ru} {
      root = r;
      cost_sc = d_ru;
      totree = std::make_shared<Tree> (root);
    };

    auto append(std::shared_ptr<PartialTree> tree, bool to_merge=false) {
      auto covered = std::make_shared<std::unordered_set<int>>();

      subtrees.push_back(std::make_pair(tree, to_merge));
      cost_sc += tree->cost_sc;
      for (auto t: tree->terms)
        terms.insert(t);
      
      if (to_merge) {
        totree->add_arc(std::make_pair(root, tree->v), tree->d_rv, trace_r, covered);

        for (auto t: tree->terms) {
          totree->add_arc(std::make_pair(tree->v, t), tree->distances_t.at(t), trace_t->at(t), covered, true, true);
        }
      }

      return covered;
    }

    auto append(std::shared_ptr<PartialTreeManager> tree, bool to_merge=false) {
      auto covered = std::make_shared<std::unordered_set<int>>();

      trees.push_back(std::make_pair(tree, to_merge));
      cost_sc += tree->cost_sc;
      for (auto t: tree->terms)
        terms.insert(t);

      // helper for subtree
      auto add_subtree = [&](std::shared_ptr<PartialTree> tree2, int u_=NONVERTEX) {
        auto r = u_==NONVERTEX? root : u_;
        auto trace_ = u_==NONVERTEX? trace_r : (*trace_u)[u_];
        totree->add_arc(std::make_pair(r, tree2->v), tree2->d_rv, trace_, covered);
        for (auto t: tree2->terms) {
          totree->add_arc(std::make_pair(tree2->v, t), tree2->distances_t.at(t), trace_t->at(t), covered, true, true);
        }
      };
      
      // helper for tree
      auto add_tree = [&](const PartialTreeManager *tree3) {
        totree->add_arc(std::make_pair(tree3->root, tree3->u), tree3->d_ru, trace_r, covered);

        for (auto tree2_flag: tree3->subtrees) {
          if (tree2_flag.second) 
            continue;
          add_subtree(tree2_flag.first, tree3->u);
        }
      };

      if (to_merge) {
        if (tree->u != NONVERTEX) {
          // could be r->{u}->{v}-{t}
          add_tree(tree.get());
        } 
        else {
          // could be r->{v}-{t}
          // appended 2-level to a manager w/ u
          for (auto tree2_flag: tree->subtrees) {
            if (tree2_flag.second)
              continue;
            add_subtree(tree2_flag.first);
          }
        }
      }

      return covered;
    }

    auto to_tree () {
      // helper for subtree
      auto add_subtree = [&](std::shared_ptr<PartialTree> tree2, int u_=NONVERTEX) {
        auto r = u_==NONVERTEX? root : u_;
        auto trace_ = u_==NONVERTEX? trace_r : (*trace_u)[u_];
        totree->add_arc(std::make_pair(r, tree2->v), tree2->d_rv, trace_);
        for (auto t: tree2->terms) {
          totree->add_arc(std::make_pair(tree2->v, t), tree2->distances_t.at(t), trace_t->at(t), nullptr, true, true);
        }
      };
      
      // helper for tree
      auto add_tree = [&](const PartialTreeManager *tree3) {
        totree->add_arc(std::make_pair(tree3->root, tree3->u), tree3->d_ru, trace_r);

        for (auto tree2_flag: tree3->subtrees) {
          if (tree2_flag.second) 
            continue;
          add_subtree(tree2_flag.first, tree3->u);
        }
      };

      if (u != NONVERTEX) { // add r->u->{v}->{t}
        add_tree(this);
      }
      else {
        for (auto tree3_flag: trees) {
          if (tree3_flag.second)
            continue;
          auto tree3 = tree3_flag.first;

          if (tree3->u != NONVERTEX) {
            // could be r->{u}->{v}-{t}
            add_tree(tree3.get());
          } 
          else {
            // could be r->{v}-{t}
            // appended 2-level to a manager w/ u
            for (auto tree2_flag: tree3->subtrees) {
              if (tree2_flag.second)
                continue;
              add_subtree(tree2_flag.first);
            }
          }
        }

        // could be r->{v}-{t}
        for (auto tree2_flag: subtrees) {
          if (tree2_flag.second)
            continue;
          add_subtree(tree2_flag.first);
        }
      }

      return totree;
    }
  };


  struct TableCell {
    int nterm;
    double d_rv_UB, cost_terms;

    TableCell(double d_rv_UB, double cost_terms, int nterm) : 
        nterm {nterm}, d_rv_UB {d_rv_UB}, cost_terms {cost_terms} {
    };
  };

  class PartialTreeTable {
    // fix v, an oracle to return terminal subsets (or their costs) given d_rv
    // delete() and partree() take O(k) time
    public:
    int v {NONVERTEX}; // the root has a single child
    std::vector<std::pair<int,double>> term_ds;
    std::vector<std::shared_ptr<TableCell>> cells; // a cell keeps a threshold 

    PartialTreeTable() {};

    PartialTreeTable(int v) : v {v} {};

    void add_term(int t, double d_vt) {
      term_ds.push_back(std::make_pair(t, d_vt));
    }

    void build(bool sorted=false) {
      if (term_ds.size() == 0)
        return;

      if (not sorted) {
        std::sort(
          term_ds.begin(),
          term_ds.end(),
          [](const std::pair<int,double> &a,  
            const std::pair<int,double> &b) { 
            return a.second < b.second;
          }
        );  
      }

      int i = 1;
      auto it = term_ds.begin(); 
      double sum = it->second, avg = it->second;
      for (it++, i++; it != term_ds.end(); it++, i++) {
        double sum_new = sum + it->second;
        double avg_new = sum_new / i;
        double thr = (avg_new - avg) * (i-1) * i; // may be zero
        cells.push_back(std::make_shared<TableCell> (thr, sum, i-1));
        sum = sum_new;
        avg = avg_new;
      }
      cells.push_back(std::make_shared<TableCell> (std::numeric_limits<double>::max(), sum, i-1));
    }

    auto find(double d_rv) {
      // fist cell that is greater than d_rv, so include as many t's as possible
      // include all terms w/ zero thr as long as d_rv > 0
      auto tar = std::make_shared<TableCell> (d_rv, 0, 0);
      auto idx = std::upper_bound(cells.begin(), cells.end(), tar,
        [](const std::shared_ptr<TableCell> a, const std::shared_ptr<TableCell> b) { 
          return a->d_rv_UB < b->d_rv_UB;
        }
      );
      if (idx == cells.end())
        return *(cells.rbegin());
      return *idx;
    }

    double density(double d_rv) {
      if (term_ds.size() == 0)
        return std::numeric_limits<double>::max();

      auto cell = find(d_rv);
      return (d_rv + cell->cost_terms) / cell->nterm;
    }

    auto partree(int root, double d_rv, const std::unordered_set<int> *terms=nullptr) {
      // no need to find(d_rv)
      auto par = std::make_shared<PartialTree> (root, v, d_rv);
      for (auto it = term_ds.begin(); it != term_ds.end(); it++) {
        if (terms != nullptr and not has_key(*terms, it->first))
          continue;
        par->add_term(it->first, it->second);
      }
      return par;
    }

    void erase(const std::unordered_set<int> &terms_del) {
      // delete v's from term_ds, and rebuild thresholds
      std::vector<std::pair<int,double>> tds;
      for (auto &p : term_ds) {
        if (has_key(terms_del, p.first))
          continue;
        tds.push_back(p);
      }

      term_ds.clear();
      term_ds = std::move(tds);
      cells.clear();
      build(true);
    }
  };


  class ThresholdedMinDensity {
    public:
    int n_thresholds;
    double thr_max;
    std::shared_ptr<std::vector<std::pair<double,double>>> thr_mindens;

    ThresholdedMinDensity(
        int n_thresholds, 
        double thr_max,
        std::unordered_map<int, PartialTreeTable> &tbls
    ) : n_thresholds {n_thresholds}, thr_max {thr_max} 
    {
      thr_mindens = std::make_shared<std::vector<std::pair<double,double>>> ();

      std::priority_queue<std::tuple<double,int,int>, 
                          std::vector<std::tuple<double,int,int>>, 
                          std::greater<std::tuple<double,int,int>>> pq_v;
      for (int i=1; i <= n_thresholds; i++) {
        double thr = thr_max/n_thresholds * i;

        // initialize pq_v with first density for each v
        if (i == 1) {
          for (auto &p: tbls) {
            auto v = p.first;
            double den = p.second.density(thr); 
            pq_v.emplace(den, v, 1);
          }
        }

        // find min density for thr
        while (true) {
          auto [den, v, i_] = pq_v.top();
          if (i == i_) {
            thr_mindens->push_back(std::make_pair(thr, den));
            break;
          } else {
            pq_v.pop();
            double den = tbls.at(v).density(thr); 
            pq_v.emplace(den, v, i);
          }
        }
      }
    };

    double min_density(double d_rv) {
      int left = 0;
      int right = thr_mindens->size() - 1;
      int idx = -1;

      // binary search the largest thr_idx <= d_rv
      while (left <= right) {
        int mid = left + (right - left) / 2;

        if (lq((*thr_mindens)[mid].first, d_rv)) {
          idx = mid;
          left = mid + 1;
        } else {
          right = mid - 1;
        }
      }

      if (idx >= 0) {
        return (*thr_mindens)[idx].second;
      }
      return 0;
    }

  };

} // namespace dst

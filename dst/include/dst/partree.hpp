#pragma once

#include <limits>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <list>

#include "dst/utils.hpp"
#include "dst/tree.hpp"


namespace dst {

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
    double distance_lastly_added = 0;
    double density_ = -1; // fake density
    bool ready = false; 

    PartialTree() {};

    PartialTree(int r, int v, double d_rv) : 
        v {v}, d_rv {d_rv} {
      root = r;
      cost_sc = d_rv;
    }

    PartialTree(double density) : density_ {density} {
      ; // construct a fake tree with a specific density
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
        distance_lastly_added = d_vt;
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
      distance_lastly_added = 0;
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
      distance_lastly_added = 0;

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

    double density_LB(int k) {
      if (terms.size() == 0)
        return -1;

      if (ready)
        return density();
      else
        return (cost_sc + (k-terms.size()) * distance_lastly_added) / k;
    }

    double density() {
      if (density_ > 0) 
        return density_;

      return BasePartialTree::density();
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
#pragma once

#include <limits>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <list>

#include <dst/tree.hpp>


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
    std::list<int> terms_after_ready;
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
        cost_sc += d_vt;
        distance_lastly_added = d_vt;
      }
    }

    void erase_and_reset(const std::unordered_set<int> &terms_del) {
      for (auto t: terms_del) {
        if (not has_key(terms, t))
          continue;
        terms.erase(t);
        cost_sc -= distances_t.at(t);
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

    auto to_tree () {
      // TODO: test
      auto tree = std::make_shared<Tree> (root);
      tree->add_arc(std::make_pair(root, v), d_rv, trace_r);
      for (auto t: terms) {
        tree->add_arc(std::make_pair(v, t), distances_t.at(t), trace_t->at(t), true, true);
      }
      return tree;
    }
  };


  class PartialTreeManager: public BasePartialTree {
    // can be used a 3-level partial tree if u != NONVERTEX
    public:
    int u {NONVERTEX}; // the root has a single child
    double d_ru;
    std::vector<std::shared_ptr<PartialTree>> subtrees; // used as r->{v}->{t}
    std::vector<std::shared_ptr<PartialTreeManager>> trees; // used as r->{u}->{v}->{t}

    PartialTreeManager() {};

    PartialTreeManager(int r) {
      root = r;
    };

    PartialTreeManager(int r, int u, double d_ru) : 
        u {u}, d_ru {d_ru} {
      root = r;
      cost_sc = d_ru;
    };

    void append(std::shared_ptr<PartialTree> tree) {
      subtrees.push_back(tree);
      cost_sc += tree->cost_sc;
      for (auto t: tree->terms)
        terms.insert(t);
    }

    void append(std::shared_ptr<PartialTreeManager> tree) {
      trees.push_back(tree);
      cost_sc += tree->cost_sc;
      for (auto t: tree->terms)
        terms.insert(t);
    }

    auto to_tree () {
      auto tree = std::make_shared<Tree> (root);

      // a helper
      auto to_tree3 = [&] (const PartialTreeManager *tree3) {
        tree->add_arc(std::make_pair(tree3->root, tree3->u), tree3->d_ru, trace_r);
        for (auto tree2: tree3->subtrees) {
          tree->add_arc(std::make_pair(tree3->u, tree2->v), tree2->d_rv, (*trace_u)[tree3->u]);
          for (auto t: tree2->terms) {
            tree->add_arc(std::make_pair(tree2->v, t), tree2->distances_t.at(t), trace_t->at(t), true, true);
          }
        }
      };

      if (u != NONVERTEX) {
        // building r->u->{v}-{t}
        to_tree3(this);
      } else {
        // building r->{u}->{v}-{t}
        for (auto tree3: trees) {
          to_tree3(tree3.get());
        }

        // building r->{v}-{t}
        for (auto tree2: subtrees) {
          tree->add_arc(std::make_pair(root, tree2->v), tree2->d_rv, trace_r);
          for (auto t: tree2->terms) {
            tree->add_arc(std::make_pair(tree2->v, t), tree2->distances_t.at(t), trace_t->at(t), true, true);
          }
        }
      }

      return tree;
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
    // delete and partree take O(k) time
    public:
    int v {NONVERTEX}; // the root has a single child
    std::vector<std::pair<int,double>> term_ds;
    std::vector<std::shared_ptr<TableCell>> cells;

    PartialTreeTable(int v) : v {v} {};

    void add_term(int t, double d_vt) {
      term_ds.push_back(std::make_pair(t, d_vt));
    }

    void build(bool sorted=false) {
      assert(term_ds.size() > 0);

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
      auto cell = find(d_rv);
      return (d_rv + cell->cost_terms) / cell->nterm;
    }

    auto partree(int root, double d_rv) {
      // no need to find(d_rv)
      auto par = std::make_shared<PartialTree> (root, v, d_rv);
      for (auto it = term_ds.begin(); it != term_ds.end(); it++) {
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

      term_ds = std::move(tds);
      cells.clear();
      build(true);
    }
  };


} // namespace dst
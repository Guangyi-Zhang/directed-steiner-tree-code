#pragma once

#include <limits>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <list>

#include <dst/tree.hpp>


namespace dst {

  class PartialTree {
    public:
    int root {NONVERTEX}; 
    int v {NONVERTEX}; // the root has a single child
    double d_rv;

    std::unordered_set<int> terms;
    std::unordered_map<int, double> distances_t;
    std::list<int> terms_after_ready;
    double cost_sc = 0;
    double distance_lastly_added = 0;
    double density_ = -1; // fake density
    bool ready = false; 

    PartialTree() {};

    PartialTree(int root, int v, double d_rv) : 
        root {root}, v {v}, d_rv {d_rv}, cost_sc {d_rv} {
      ;
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
        if (ready) 
          terms_after_ready.push_front(t);
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

      if (terms.size() == 0)
        return std::numeric_limits<double>::max();
      return cost_sc / terms.size();
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

    auto to_tree (
        const std::shared_ptr<std::unordered_map<int,int>> trace_r, 
        const std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>> trace_t
    ) {
      auto tree = std::make_shared<Tree> (root);
      tree->add_arc(std::make_pair(root, v), d_rv, trace_r);
      for (auto t: terms) {
        tree->add_arc(std::make_pair(v, t), distances_t.at(t), trace_t->at(t), true, true);
      }
      return tree;
    }
  };


  class PartialTreeManager {
    // can be used a 3-level partial tree if u != NONVERTEX
    public:
    int root {NONVERTEX}; 
    int u {NONVERTEX}; // the root has a single child
    double d_ru;
    std::vector<std::shared_ptr<PartialTree>> subtrees;

    double cost_sc = 0;
    std::unordered_set<int> terms;

    PartialTreeManager() {};

    PartialTreeManager(int r) : root {r} {};

    PartialTreeManager(int r, int u, double d_ru) : 
        root {r}, u {u}, d_ru {d_ru}, cost_sc {d_ru} {};

    double density() {
      if (terms.size() == 0)
        return std::numeric_limits<double>::max();
      return cost_sc / terms.size();
    }

    void append(std::shared_ptr<PartialTree> tree) {
      subtrees.push_back(tree);
      cost_sc += tree->cost_sc;
      for (auto t: tree->terms)
        terms.insert(t);
    }

    auto to_tree (
        const std::shared_ptr<std::unordered_map<int,int>> trace_r, 
        const std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>> trace_t,
        const std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>> trace_u = nullptr
    ) {
      auto tree = std::make_shared<Tree> (root);

      if (u != NONVERTEX) {
        tree->add_arc(std::make_pair(root, u), d_ru, trace_r);
        // TODO
      } else {
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

} // namespace dst
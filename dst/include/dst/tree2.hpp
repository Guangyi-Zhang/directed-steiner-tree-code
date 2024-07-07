#pragma once

#include <limits>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <list>

#include <dst/tree.hpp>


namespace dst {

  class Level2PartialTree {
    public:
    int root {NONVERTEX}; 
    int v {NONVERTEX}; // the root has a single child
    double d_rv;

    std::unordered_set<int> terms;
    std::unordered_map<int, double> distances_v; // b/w v and terms
    std::list<int> terms_after_ready;
    double cost_sc = 0;
    double distance_lastly_added = 0;
    double density_ = -1; // fake density
    bool ready = false; 

    Level2PartialTree() {};

    Level2PartialTree(int root, int v, double d_rv) : 
        root {root}, v {v}, d_rv {d_rv}, cost_sc {d_rv} {
      ;
    }

    Level2PartialTree(double density) : density_ {density} {
      ; // construct a fake tree with a specific density
    }

    void add_term(int t, double d_vt) {
      distances_v[t] = d_vt;
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
        cost_sc -= distances_v.at(t);
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

        add_term(t, distances_v.at(t));
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
      auto tree = std::make_shared<PartialTree> (root);
      tree->add_arc(std::make_pair(root, v), d_rv, trace_r);
      for (auto t: terms) {
        tree->add_arc(std::make_pair(v, t), distances_v.at(t), trace_t->at(t), true, true);
      }
      return tree;
    }

  };

} // namespace dst
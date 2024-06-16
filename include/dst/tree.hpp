#pragma once

#include <limits>
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>

#include <dst/consts.hpp>
#include <dst/utils.hpp>


namespace dst {

  class PartialTree {
  public:
    static const std::map<std::pair<int,int>, double>* w; 

    int root;
    double cost {0}, cost_sc {0};
    std::unordered_map<int,int> trace_sc; // trace "short-cut" chosen 2-level trees
    std::unordered_map<int,int> trace;
    std::unordered_set<int> terms_cov;

    PartialTree(int root) : 
        root {root} {
      trace[root] = NONVERTEX;
      trace_sc[root] = NONVERTEX;
    }

    PartialTree(std::pair<int,int> arc, double w_arc, 
        const std::unordered_map<int,int> &trace_arc, 
        bool reverse=false,
        bool is_terminal=false) : 
        root {arc.first} {
      if (is_terminal)
        terms_cov.insert(arc.second);

      trace[arc.first] = NONVERTEX;
      trace_sc[arc.first] = NONVERTEX;
      trace[arc.second] = arc.first;
      trace_sc[arc.second] = arc.first;

      cost_sc += w_arc;

      auto v = reverse? arc.first: arc.second;
      while (trace_arc.at(v) != NONVERTEX) {
        auto e = reverse? std::make_pair(v, trace_arc.at(v)) : 
                          std::make_pair(trace_arc.at(v), v);
        cost += PartialTree::w->at(e);
        trace[e.second] = e.first;
        v = trace_arc.at(v);
      }
    }

    void append(PartialTree &tree) {
      // append either r-u-{v}-{t}, or u-{v}-{r}
      cost_sc += tree.cost_sc; // ok to count (u,v) multi times
      auto es = edges();
      for (auto &e: tree.edges()) {
        if (not has_key(es, e)) {
          cost += PartialTree::w->at(e);
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
        //if (has_key(trace, p.first))
        //  assert (p.second == trace[p.first]);

        // TODO: i think a bug here
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

  const std::map<std::pair<int,int>, double>* PartialTree::w = nullptr; 
}
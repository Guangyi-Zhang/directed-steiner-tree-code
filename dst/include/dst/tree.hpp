#pragma once

#include <limits>
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>

#include <fmt/ranges.h>
#include <boost/functional/hash.hpp>
#include <dst/consts.hpp>
#include <dst/utils.hpp>


namespace dst {


  class TreeNode {
  public:
    int v {NONVERTEX}, pa {NONVERTEX};
    std::unordered_set<int> children;

    TreeNode() {
      ;
    }

    TreeNode(int v) : v {v} {
      ;
    }
  };


  class PartialTree {
  public:
    static const std::unordered_map<std::pair<int,int>, double, boost::hash<std::pair<int,int>>>* w; 

    int root;
    double cost {0}, cost_sc {0};
    std::unordered_map<int,int> trace_sc; // trace "short-cut" chosen 2-level trees
    std::unordered_map<int,int> trace;
    std::unordered_set<int> terms_cov;

    std::unordered_map<std::string, std::string> debuginfo;

    PartialTree(int root) : 
        root {root} {
      trace[root] = NONVERTEX;
      trace_sc[root] = NONVERTEX;
    }

    void add_arc(std::pair<int,int> arc, double w_arc, 
        const std::unordered_map<int,int> &trace_arc, 
        bool reverse=false,
        bool is_terminal=false) {
      // start from PartialTree(root), arcs must be added from top-down
      // i.e., arc.first must exist already beforehand
      if (is_terminal)
        terms_cov.insert(arc.second);

      trace_sc[arc.second] = arc.first;
      cost_sc += w_arc;

      auto v = reverse? arc.first: arc.second;
      while (trace_arc.at(v) != NONVERTEX) {
        auto e = reverse? std::make_pair(v, trace_arc.at(v)) : 
                          std::make_pair(trace_arc.at(v), v);
        if (not has_key(trace, e.second)) {
          cost += PartialTree::w->at(e);
          trace[e.second] = e.first;
        }
        else if (reverse) {
          ; // need to prune disconnected dangling branches later
        }
        else {
          // this is fine to have
          // e.first == NONVERTEX or e.first == trace[e.second].
          // there exist cycles containing u->root->v or non-tree crosing,
          // ok to shortcut 
          break;
        }
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
        if (not has_key(trace, p.first))
          trace[p.first] = p.second;
        else {
          if (p.second == NONVERTEX or p.second == trace[p.first]) {
            // this is fine
          } else {
            // there exist cycles containing u->root->v or non-tree crosing,
            // ok to shortcut but remember to prune disconnected dangling arcs w/o terminals
            if (DEBUG) fmt::println("append: arc ({},{}) while arc ({},{}) in tree", p.second, p.first, trace[p.first], p.first);
          }
        }
      }
    }

    double density() const {
      if (terms_cov.size() == 0)
        return std::numeric_limits<double>::max();
      return cost_sc / terms_cov.size();
    }

    bool zero_coverage() const {
      return terms_cov.size() == 0 ? true: false;
    }

    std::unordered_set<std::pair<int,int>, boost::hash<std::pair<int,int>>> edges() const {
      std::unordered_set<std::pair<int,int>, boost::hash<std::pair<int,int>>> es;

      for (auto t: terms_cov) {
        int v {t};
        while (trace.at(v) != NONVERTEX) {
          es.insert(std::make_pair(trace.at(v), v));
          v = trace.at(v);
        }
      }

      return es;
    }

    std::unordered_map<int, TreeNode> to_treenode() const {
      std::unordered_map<int, TreeNode> v2nd;

      for (auto t: terms_cov) {
        int v {t};
        while (trace.at(v) != NONVERTEX) {
          int pv = trace.at(v);
          if (not has_key(v2nd, pv))
            v2nd[pv] = TreeNode {pv};
          v2nd.at(pv).children.insert(v);
          if (not has_key(v2nd, v))
            v2nd[v] = TreeNode {v};
          v2nd.at(v).pa = v;

          v = trace.at(v);
        }
      }

      return v2nd;
    }

    void print() {
      auto &&v2nd = to_treenode();
      if (v2nd.empty())
        return;
      _print_treenode(root, NONVERTEX, v2nd, "", true);
    }

    void _print_treenode(int v, int pa, std::unordered_map<int, TreeNode> &v2nd, 
        std::string prefix, bool is_last_child) {
      auto &nd = v2nd.at(v);
      fmt::print(prefix);
      fmt::print(is_last_child ? "└──" : "├──");
      if (pa != NONVERTEX)
        fmt::println("{} ({:05.1f})", v, PartialTree::w->at(std::make_pair(pa, v)));
      else
        fmt::println("{}", v);
      size_t i = 0;
      for (auto child : nd.children) {
      //for (size_t i=0; i < nd.children.size(); i++) {
        _print_treenode(child, v, v2nd, 
            prefix + (is_last_child ? "    " : "│   "),
            i == nd.children.size()-1);
        i++;
      }
    }

    double cost_trimmed() const {
      // exclude dangling arcs
      double sum {0};
      auto &&es = edges();
      for (const auto &e: es) {
        sum += PartialTree::w->at(e);
      }
      return sum;
    }
  };

  const std::unordered_map<std::pair<int,int>, double, boost::hash<std::pair<int,int>>>* PartialTree::w = nullptr; 
}

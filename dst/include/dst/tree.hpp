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
#include "dst/consts.hpp"
#include "dst/utils.hpp"
#include "dst/csrgraph.hpp"


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


  class Tree {
  public:
    static const CSRGraph* csr_graph;

    int root;
    double cost {0}, cost_sc {0};
    std::unordered_map<int,int> trace;
    std::unordered_set<int> terms_cov;

    std::unordered_map<std::string, std::string> debuginfo;

    Tree() {};

    Tree(int root) : 
        root {root} {
      trace[root] = NONVERTEX;
    }

    void add_arc(
        std::pair<int,int> arc, 
        double w_arc, 
        const std::shared_ptr<std::unordered_map<int,int>> trace_arc, 
        std::shared_ptr<std::unordered_set<int>> covered=nullptr,
        bool reverse=false,
        bool is_terminal=false) {
      // start from Tree(root), arcs must be added from top-down
      // i.e., arc.first must exist already beforehand
      if (arc.first == arc.second)
        return;

      if (is_terminal)
        terms_cov.insert(arc.second);

      cost_sc += w_arc;

      auto v = reverse? arc.first: arc.second;
      while (trace_arc->at(v) != NONVERTEX) {
        auto e = reverse? std::make_pair(v, trace_arc->at(v)) : 
                          std::make_pair(trace_arc->at(v), v);
        if(DEBUG) fmt::println("{}<-{}, {}", e.first, e.second, not has_key(trace, e.second));
        if (not has_key(trace, e.second)) {
          cost += Tree::csr_graph->edge_weight(e.first, e.second);
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
        if (covered != nullptr) 
          covered->insert(v);
        v = trace_arc->at(v);
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

    std::shared_ptr<std::unordered_set<std::pair<int,int>, boost::hash<std::pair<int,int>>>> 
    edges() const
    {
      auto es = std::make_shared<std::unordered_set<std::pair<int,int>, boost::hash<std::pair<int,int>>>> ();

      for (auto t: terms_cov) {
        int v {t};
        while (trace.at(v) != NONVERTEX) {
          es->insert(std::make_pair(trace.at(v), v));
          v = trace.at(v);
        }
      }

      return es;
    }

    std::unordered_map<int, TreeNode> to_treenode() const 
    {
      // TODO: made into pointer
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

    // Verify that each required terminal is covered by a root-directed path.
    bool is_valid(const std::unordered_set<int> &required_terms) const {
      if (Tree::csr_graph == nullptr ||
          root < 0 || root >= Tree::csr_graph->n ||
          !has_key(trace, root) || trace.at(root) != NONVERTEX)
        return false;

      for (const auto terminal : required_terms) {
        if (!has_key(terms_cov, terminal) || !has_key(trace, terminal))
          return false;

        std::unordered_set<int> path_vertices;
        int v = terminal;
        while (v != root) {
          if (!path_vertices.insert(v).second || !has_key(trace, v))
            return false;  // Cycle or a path that stops before the root.

          const int parent = trace.at(v);
          if (parent < 0 || parent >= Tree::csr_graph->n ||
              v < 0 || v >= Tree::csr_graph->n ||
              Tree::csr_graph->edge_weight(parent, v) == std::numeric_limits<float>::max())
            return false;
          v = parent;
        }
      }

      return true;
    }

    void _print_treenode(
        int v, 
        int pa, 
        std::unordered_map<int, TreeNode> &v2nd, 
        std::string prefix, 
        bool is_last_child
    ) {
      auto &nd = v2nd.at(v);
      fmt::print(prefix);
      fmt::print(is_last_child ? "└──" : "├──");
      if (pa != NONVERTEX)
        fmt::println("{} ({:05.1f})", v, Tree::csr_graph->edge_weight(pa, v));
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
      auto es = edges();
      for (const auto &e: *es) {
        sum += Tree::csr_graph->edge_weight(e.first, e.second);
      }
      return sum;
    }
  };

  const CSRGraph* Tree::csr_graph = nullptr;

} // namespace dst

#pragma once

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <queue>
#include <vector>
#include <list>
#include <functional> // std::greater

#include <boost/container_hash/hash.hpp>
#include <dst/consts.hpp>


namespace dst {

  auto dijkstra(
      const std::unordered_map<int, std::vector<int>> &adj, 
      const std::unordered_map<std::pair<int,int>, double, boost::hash<std::pair<int,int>>> &edgeweight, 
      int source,
      bool reverse=false,
      const std::unordered_set<int> &to_reach=std::unordered_set<int> {}
  ) {
    auto distances = std::make_shared<std::unordered_map<int,double>>();
    auto trace = std::make_shared<std::unordered_map<int,int>>();
    std::priority_queue<std::tuple<double,int,int>, 
                        std::vector<std::tuple<double,int,int>>, 
                        std::greater<std::tuple<double,int,int>>> pq;
    pq.emplace(0, NONVERTEX, source);
    size_t n_reached = 0;

    while (!pq.empty()) {
      double d_u;
      int u_prev, u; 
      std::tie(d_u, u_prev, u) = pq.top();
      pq.pop();

      if (has_key(*distances, u)) 
        continue;
      (*distances)[u] = d_u;
      (*trace)[u] = u_prev;
      if (to_reach.size() > 0 and has_key(to_reach, u)) {
        n_reached++;
        if (to_reach.size() == n_reached)
          break;
      }

      if (not has_key(adj, u))
        continue;
      for (const auto& v: adj.at(u)) {
        if (has_key(*distances, v)) 
          continue;
        double weight = reverse? edgeweight.at({v,u}) : edgeweight.at({u,v});
        pq.emplace(d_u + weight, u, v);
      }
    }

    return std::make_tuple(distances, trace);
  }


  class Dijkstra {
    public:
    const std::unordered_map<int, std::vector<int>> *padj {nullptr};
    const std::unordered_map<std::pair<int,int>, double, boost::hash<std::pair<int,int>>> *pedgeweight {nullptr};
    int source;
    bool reverse=false;

    std::shared_ptr<std::unordered_map<int,double>> distances;
    std::shared_ptr<std::unordered_map<int,int>> trace;
    std::priority_queue<std::tuple<double,int,int>, 
                        std::vector<std::tuple<double,int,int>>, 
                        std::greater<std::tuple<double,int,int>>> pq;

    Dijkstra() {};

    Dijkstra(
        const std::unordered_map<int, std::vector<int>> &adj, 
        const std::unordered_map<std::pair<int,int>, double, boost::hash<std::pair<int,int>>> &edgeweight, 
        int source,
        bool reverse=false) :
        padj {&adj}, pedgeweight {&edgeweight}, source {source}, reverse {reverse} 
    {
      distances = std::make_shared<std::unordered_map<int,double>>();
      trace = std::make_shared<std::unordered_map<int,int>>();
      pq.emplace(0, NONVERTEX, source);
    }

    std::tuple<int,double> next() {
      while (not pq.empty()) {
        double d_u;
        int u_prev, u; 
        std::tie(d_u, u_prev, u) = pq.top();
        pq.pop();

        if (has_key(*distances, u)) 
          continue;

        (*distances)[u] = d_u;
        (*trace)[u] = u_prev;

        if (not has_key(*padj, u))
          continue;

        for (const auto& v: padj->at(u)) {
          if (has_key(*distances, v)) 
            continue;
          double weight = reverse? pedgeweight->at({v,u}) : pedgeweight->at({u,v});
          pq.emplace(d_u + weight, u, v);
        }

        return std::make_tuple(u, d_u); // break
      }

      // nothing to return
      return std::make_tuple(NONVERTEX, -1);
    }

  }; // end of class


  /*
    Coordindated SSSP's, one from each terminal
  */
  class CoordinatedDijkstra {
    public:
    const std::unordered_map<int, std::vector<int>> *padj {nullptr};
    const std::unordered_map<std::pair<int,int>, double, boost::hash<std::pair<int,int>>> *pedgeweight {nullptr};
    const std::unordered_set<int> *psources {nullptr};
    bool reverse=false;

    std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,double>>>> distances_t;
    std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>> trace_t;
    std::unordered_set<int> sources_del;
    std::priority_queue<std::tuple<double,int,int,int>, 
                        std::vector<std::tuple<double,int,int,int>>, 
                        std::greater<std::tuple<double,int,int,int>>> pq;


    CoordinatedDijkstra(
        const std::unordered_map<int, std::vector<int>> &adj, 
        const std::unordered_map<std::pair<int,int>, double, boost::hash<std::pair<int,int>>> &edgeweight, 
        const std::unordered_set<int> &sources,
        bool reverse=false) :
        padj {&adj}, pedgeweight {&edgeweight}, psources {&sources}, reverse {reverse} 
    {
      distances_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,double>>>>();
      trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      for (auto &source: *psources) {
        pq.emplace(0, NONVERTEX, source, source);
      }
    }


    void delete_source(int src) {
      sources_del.insert(src);
    }


    std::tuple<int,int,double> next() {
      while (not pq.empty()) {
        double d_u;
        int u_prev, u, source; 
        std::tie(d_u, u_prev, u, source) = pq.top();
        pq.pop();

        if (has_key(sources_del, source)) 
          continue;
        if (not has_key(*distances_t, source)) {
          (*distances_t)[source] = std::make_shared<std::unordered_map<int,double>>();
          (*trace_t)[source] = std::make_shared<std::unordered_map<int,int>>();
        }
        auto &distances = distances_t->at(source);
        if (has_key(*distances, u)) 
          continue;

        (*distances)[u] = d_u;
        (*(trace_t->at(source)))[u] = u_prev;
        if (has_key(*padj, u)) {
          for (const auto& v: padj->at(u)) {
            if (has_key(*distances, v)) 
              continue;
            double weight = reverse? pedgeweight->at({v,u}) : pedgeweight->at({u,v});
            pq.emplace(d_u + weight, u, v, source);
          }
        }

        return std::make_tuple(source, u, d_u); // break
      }

      // nothing to return
      return std::make_tuple(NONVERTEX, NONVERTEX, -1);
    }

  }; // end of class

}
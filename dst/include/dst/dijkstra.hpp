#pragma once

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <queue>
#include <vector>
#include <list>
#include <functional> // std::greater

#include "dst/consts.hpp"


namespace dst {

  auto dijkstra(
      const std::vector<int> &offset,
      const std::vector<int> &target,
      const std::vector<float> &weight,
      int source,
      const std::unordered_set<int> &to_reach=std::unordered_set<int> {}
  ) {
    auto distances = std::make_shared<std::unordered_map<int,float>>();
    auto trace = std::make_shared<std::unordered_map<int,int>>();
    std::priority_queue<std::tuple<float,int,int>,
                        std::vector<std::tuple<float,int,int>>,
                        std::greater<std::tuple<float,int,int>>> pq;
    pq.emplace(0.0f, NONVERTEX, source);
    size_t n_reached = 0;

    while (!pq.empty()) {
      float d_u;
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

      for (int i = offset[u]; i < offset[u+1]; i++) {
        int v = target[i];
        if (has_key(*distances, v))
          continue;
        pq.emplace(d_u + weight[i], u, v);
      }
    }

    return std::make_tuple(distances, trace);
  }


  // Shortest paths from `sources` (dist=0 at each source).
  auto dijkstra_multisource(
      const std::vector<int> &offset,
      const std::vector<int> &target,
      const std::vector<float> &weight,
      const std::unordered_set<int> &sources,
      const std::unordered_set<int> &to_reach=std::unordered_set<int> {}
  ) {
    auto distances = std::make_shared<std::unordered_map<int,float>>();
    auto trace = std::make_shared<std::unordered_map<int,int>>();
    // (distance, predecessor, vertex); min-heap by distance
    std::priority_queue<std::tuple<float,int,int>,
                        std::vector<std::tuple<float,int,int>>,
                        std::greater<std::tuple<float,int,int>>> pq;
    for (auto source: sources) {
      pq.emplace(0.0f, NONVERTEX, source);
    }
    size_t n_reached = 0;

    while (!pq.empty()) {
      float d_u;
      int u_prev, u;
      std::tie(d_u, u_prev, u) = pq.top();
      pq.pop();

      if (has_key(*distances, u))
        continue; // already settled with a shorter path
      (*distances)[u] = d_u;
      (*trace)[u] = u_prev;
      if (to_reach.size() > 0 and has_key(to_reach, u)) {
        n_reached++;
        if (to_reach.size() == n_reached)
          break; // stop once every target is settled
      }

      for (int i = offset[u]; i < offset[u+1]; i++) {
        int v = target[i];
        if (has_key(*distances, v))
          continue;
        pq.emplace(d_u + weight[i], u, v);
      }
    }

    return std::make_tuple(distances, trace);
  }


  class Dijkstra {
    public:
    const std::vector<int> *poffset {nullptr};
    const std::vector<int> *ptarget {nullptr};
    const std::vector<float> *pweight {nullptr};
    int source;

    std::shared_ptr<std::unordered_map<int,float>> distances;
    std::shared_ptr<std::unordered_map<int,int>> trace;
    std::priority_queue<std::tuple<float,int,int>,
                        std::vector<std::tuple<float,int,int>>,
                        std::greater<std::tuple<float,int,int>>> pq;

    Dijkstra() {};

    Dijkstra(
        const std::vector<int> &offset,
        const std::vector<int> &target,
        const std::vector<float> &weight,
        int source) :
        poffset {&offset}, ptarget {&target}, pweight {&weight}, source {source}
    {
      distances = std::make_shared<std::unordered_map<int,float>>();
      trace = std::make_shared<std::unordered_map<int,int>>();
      pq.emplace(0.0f, NONVERTEX, source);
    }

    std::tuple<int,float> next() {
      while (not pq.empty()) {
        float d_u;
        int u_prev, u;
        std::tie(d_u, u_prev, u) = pq.top();
        pq.pop();

        if (has_key(*distances, u))
          continue;

        (*distances)[u] = d_u;
        (*trace)[u] = u_prev;

        for (int i = (*poffset)[u]; i < (*poffset)[u+1]; i++) {
          int v = (*ptarget)[i];
          if (has_key(*distances, v))
            continue;
          pq.emplace(d_u + (*pweight)[i], u, v);
        }

        return std::make_tuple(u, d_u);
      }

      return std::make_tuple(NONVERTEX, -1.0f);
    }

  }; // end of class


  /*
    Coordindated SSSP's, one from each terminal
  */
  class CoordinatedDijkstra {
    public:
    const std::vector<int> *poffset {nullptr};
    const std::vector<int> *ptarget {nullptr};
    const std::vector<float> *pweight {nullptr};
    const std::unordered_set<int> *psources {nullptr};

    std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,float>>>> distances_t;
    std::shared_ptr<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>> trace_t;
    std::unordered_set<int> sources_del;
    std::priority_queue<std::tuple<float,int,int,int>,
                        std::vector<std::tuple<float,int,int,int>>,
                        std::greater<std::tuple<float,int,int,int>>> pq;


    CoordinatedDijkstra(
        const std::vector<int> &offset,
        const std::vector<int> &target,
        const std::vector<float> &weight,
        const std::unordered_set<int> &sources) :
        poffset {&offset}, ptarget {&target}, pweight {&weight}, psources {&sources}
    {
      distances_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,float>>>>();
      trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
      for (auto &source: *psources) {
        pq.emplace(0.0f, NONVERTEX, source, source);
      }
    }


    void delete_source(int src) {
      sources_del.insert(src);
    }


    std::tuple<int,int,float> next() {
      while (not pq.empty()) {
        float d_u;
        int u_prev, u, source;
        std::tie(d_u, u_prev, u, source) = pq.top();
        pq.pop();

        if (has_key(sources_del, source))
          continue;
        if (not has_key(*distances_t, source)) {
          (*distances_t)[source] = std::make_shared<std::unordered_map<int,float>>();
          (*trace_t)[source] = std::make_shared<std::unordered_map<int,int>>();
        }
        auto &distances = distances_t->at(source);
        if (has_key(*distances, u))
          continue;

        (*distances)[u] = d_u;
        (*(trace_t->at(source)))[u] = u_prev;
        for (int i = (*poffset)[u]; i < (*poffset)[u+1]; i++) {
          int v = (*ptarget)[i];
          if (has_key(*distances, v))
            continue;
          pq.emplace(d_u + (*pweight)[i], u, v, source);
        }

        return std::make_tuple(source, u, d_u);
      }

      // nothing to return
      return std::make_tuple(NONVERTEX, NONVERTEX, -1.0f);
    }

  }; // end of class

}

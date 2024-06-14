#include <doctest/doctest.h>
#include <fmt/ranges.h>
#include "dst/dst.hpp"
#include "dst/version.hpp"

#include <string>
#include <map>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <tuple>



TEST_CASE("level2_alg") {
  using namespace dst;

  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(0,2), 
                                         std::make_pair(0,3), 
                                         std::make_pair(0,4), 
                                         std::make_pair(1,11), 
                                         std::make_pair(2,12), 
                                         std::make_pair(3,13), 
                                         std::make_pair(4,11), 
                                         std::make_pair(4,12), 
                                         std::make_pair(4,13)};
  std::vector<double> weights {1,1,1,2, 2,2,2, 1,1,1};
  std::vector<int> terms {11,12,13};
  DST dt = DST(edges, weights, 0, terms);

  CHECK(dt.naive_alg() == 1*3 + 2*3);
  auto tree = dt.level2_alg();
  CHECK(tree.cost == 2 + 1*3);
  CHECK(tree.cost_sc == 2 + 1*3);
  CHECK((tree.terms_cov == std::unordered_set<int> {14,15,16}));
}


TEST_CASE("DST") {
  using namespace dst;

  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(0,2), 
                                         std::make_pair(0,3), 
                                         std::make_pair(1,2), 
                                         std::make_pair(3,4), 
                                         std::make_pair(4,3), 
                                         std::make_pair(3,5)};
  std::vector<double> weights {1,1,1,1,1,1,1};
  std::vector<int> terms {2,3,5};
  DST dt = DST(edges, weights, 0, terms);

  CHECK(dt.terms_dm == (std::unordered_set<int> {6,7,8}));
  
  double cost = dt.naive_alg();
  CHECK(cost == 1+1+2-1);
}


TEST_CASE("dijkstra") {
  using namespace dst;

  /*
    +----0----+
    |    |    |
    |    |    |
    |    |    |
    1----2    3
              |
              |
              |
              4
  */
  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(0,2), 
                                         std::make_pair(0,3), 
                                         std::make_pair(1,2), 
                                         std::make_pair(3,4), 
                                         std::make_pair(4,3)};
  std::vector<double> weights {1,1,1,1,1,1};
  std::vector<int> terms {};
  DST dt = DST(edges, weights, 0, terms);

  auto p1 = dijkstra(dt.adj, dt.w, 0);
  auto dists = p1.first;
  auto trace = p1.second;
  CHECK(dists[0] == 0);
  CHECK(dists[1] == 1);
  CHECK(dists[2] == 1);
  CHECK(dists[3] == 1);
  CHECK(dists[4] == 2);
  CHECK(trace[4] == 3);

  auto p2 = dijkstra(dt.adj, dt.w, 4);
  dists = p2.first;
  trace = p2.second;
  CHECK (dists.find(0) == dists.end());
  CHECK (dists.find(1) == dists.end());
  CHECK (dists.find(2) == dists.end());
  CHECK (dists[3] == 1);
  CHECK (dists[4] == 0);

  auto p3 = dijkstra(dt.adj_r, dt.w, 4, true);
  dists = p3.first;
  trace = p3.second;
  CHECK (dists[0] == 2);
  CHECK (dists.find(1) == dists.end());
  CHECK (dists.find(2) == dists.end());
  CHECK (dists[3] == 1);
  CHECK (dists[4] == 0);
  CHECK (trace[4] == NONVERTEX);
  CHECK (trace[3] == 4);
  CHECK (trace.find(2) == trace.end());
}


TEST_CASE("misc") {
  using namespace dst;

  std::unordered_set<int> st {1,2,3};
  std::map<int,std::string> m {std::make_pair(1,"aa"), 
                               std::make_pair(2,"bb")};
  CHECK(has_key(st, 1));
  CHECK(not has_key(st, 4));
  CHECK(has_key(m, 2));
  CHECK(not has_key(m, 3));
}

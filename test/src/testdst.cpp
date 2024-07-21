#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <tuple>

#include <doctest/doctest.h>
#include <fmt/ranges.h>
#include "dst/dst.hpp"
#include "dst/dijkstra.hpp"
#include "dst/utils.hpp"

#include "tests/test_level2.cpp"
#include "tests/test_level3.cpp"


TEST_CASE("adaptive_level1") {
  using namespace dst;

  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(0,2), 
                                         std::make_pair(1,2)};
  std::vector<double> weights {1,1,0.9};
  std::vector<int> terms {1,2};
  DST dt = DST(edges, weights, 0, terms);

  CHECK(dt.level1_alg()->cost == 1+1);
  CHECK(dt.level1_alg()->cost_sc == 1+1);

  auto tree = dt.adaptive_level1_alg();
  CHECK(tree->cost == 1.9);
  CHECK(tree->cost_sc == 1.9);
}


TEST_CASE("level1") {
  /*
    <---+0+--->
    |    +    |
    |    |    |
    v    v    v
    1+-->2    3 -> 5
              +
              |
              +
              4
  */
  using namespace dst;

  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(0,2), 
                                         std::make_pair(0,3), 
                                         std::make_pair(1,2), 
                                         std::make_pair(3,4), 
                                         std::make_pair(4,3), 
                                         std::make_pair(3,5)};
  std::vector<double> weights {1,1,1,0.9,1,1,1};
  std::vector<int> terms {2,3,5};
  DST dt = DST(edges, weights, 0, terms);

  CHECK(dt.terms_dm == (std::unordered_set<int> {6,7,8}));
  
  CHECK(dt.level1_alg()->cost == 1+1+2-1);
  CHECK(dt.level1_alg()->cost_sc == 1+1+2);

  auto tree = dt.adaptive_level1_alg();
  CHECK(tree->cost == 1+1+2-1);
  CHECK(tree->cost_sc == 1+1+1);
}


TEST_CASE("CoordinatedDijkstra") {
  using namespace dst;
  /*
    <-----+0+---->
    |            |
    v            v
    1            2+-------->
    +            +         |
    |            |         |
    v            v         v
    11+--------->12        13
  */

  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(0,2), 
                                         std::make_pair(1,11), 
                                         std::make_pair(2,12), 
                                         std::make_pair(2,13), 
                                         std::make_pair(11,12)};
  std::vector<double> weights {1,1, 2,4,2.5, 1};
  std::vector<int> terms {11,12,13};
  DST dt = DST(edges, weights, 0, terms);

  std::unordered_set<int> sources {terms.begin(), terms.end()};
  CoordinatedDijkstra cosssp {dt.adj_r, dt.w, sources, true};
  int source, u; 
  double d_u;
  std::tie(source, u, d_u) = cosssp.next();
  while(eq(d_u, 0))
    std::tie(source, u, d_u) = cosssp.next();

  CHECK ((source == 12 and u == 11 and eq(d_u,1)));

  std::tie(source, u, d_u) = cosssp.next();
  CHECK ((source == 11 and u == 1 and eq(d_u,2)));

  std::tie(source, u, d_u) = cosssp.next();
  CHECK ((source == 13 and u == 2 and eq(d_u,2.5)));

  cosssp.delete_source(11);
  std::tie(source, u, d_u) = cosssp.next();
  CHECK ((source == 12 and u == 1 and eq(d_u,3)));

  std::tie(source, u, d_u) = cosssp.next();
  CHECK ((source == 13 and u == 0 and eq(d_u,3.5)));

  std::tie(source, u, d_u) = cosssp.next();
  CHECK ((source == 12 and u == 0 and eq(d_u,4)));
}



TEST_CASE("dijkstra") {
  using namespace dst;

  /*
    <---+0+--->
    |    +    |
    |    |    |
    v    v    v
    1+-->2    3
              +
              |
              +
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

  auto [dists, trace] = dijkstra(dt.adj, dt.w, 0);
  CHECK(dists->at(0) == 0);
  CHECK(dists->at(1) == 1);
  CHECK(dists->at(2) == 1);
  CHECK(dists->at(3) == 1);
  CHECK(dists->at(4) == 2);
  CHECK(trace->at(4) == 3);

  std::tie(dists, trace) = dijkstra(dt.adj, dt.w, 4);
  CHECK (dists->find(0) == dists->end());
  CHECK (dists->find(1) == dists->end());
  CHECK (dists->find(2) == dists->end());
  CHECK (dists->at(3) == 1);
  CHECK (dists->at(4) == 0);

  std::tie(dists, trace) = dijkstra(dt.adj_r, dt.w, 4, true);
  CHECK (dists->at(0) == 2);
  CHECK (dists->find(1) == dists->end());
  CHECK (dists->find(2) == dists->end());
  CHECK (dists->at(3) == 1);
  CHECK (dists->at(4) == 0);
  CHECK (trace->at(4) == NONVERTEX);
  CHECK (trace->at(3) == 4);
  CHECK (trace->find(2) == trace->end());

  auto p4 = dijkstra(dt.adj, dt.w, 0, false, {3});
  CHECK (not has_key(*std::get<0>(p4), 4));
}


TEST_CASE("pq") {
  using namespace dst;
  
  std::priority_queue<int> pq;
  std::vector<int> &tasks = Container(pq);

  for(int i=0;i<20;i++){
    pq.push(i);
  }
  CHECK(tasks.size() == pq.size());

  int sum = 0;
  for(auto i=tasks.begin(); i!=tasks.end(); i++)
    sum += *i;
  CHECK(sum == (0+19)*10);
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

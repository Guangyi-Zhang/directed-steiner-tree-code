#include <doctest/doctest.h>
#include <fmt/ranges.h>
#include "dst/dst.hpp"
#include "dst/dijkstra.hpp"
#include "dst/utils.hpp"

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <tuple>


TEST_CASE("Level2PartialTree") {
  using namespace dst;

  int root {0}, v {1};
  double d_rv {1};
  int k = 5;
  Level2PartialTree tree {root, v, d_rv};
  CHECK(not tree.is_ready());
  CHECK(tree.density() > 1e8);
  CHECK(tree.density_LB(k) < 0);

  tree.add_term(11, 1); // density = 2/1
  CHECK(eq(tree.density_LB(k), (1.0+k)/k));
  tree.add_term(12, 2); // density = 4/2
  CHECK(not tree.is_ready());
  CHECK(eq(tree.density_LB(k), (1.0+1+(k-1)*2)/k));
  tree.add_term(13, 3); // density = 7/3
  CHECK(tree.is_ready());
  CHECK(eq(tree.density(), 2));
  CHECK(eq(tree.density_LB(k), tree.density()));
  CHECK(tree.terms.size() == 2);

  tree.erase_and_reset({11, 12});
  CHECK(not tree.is_ready());
  CHECK(eq(tree.density(), 1+3));
  CHECK(eq(tree.density_LB(k), (1.0+k*3)/k));
  CHECK(tree.terms.size() == 1);

  tree.add_term(14, 4); // density = (1+3+4)/2
  CHECK(eq(tree.density_LB(k), (1.0+3+(k-1)*4)/k));
  tree.add_term(14, 5); // density = (1+3+4+5)/3
  CHECK(tree.is_ready());
  CHECK(eq(tree.density(), 4));
  CHECK(eq(tree.density_LB(k), tree.density()));
  CHECK(tree.terms.size() == 2);

  Level2PartialTree tree2 {0.5};
  CHECK(eq(tree2.density(), 0.5));
}


TEST_CASE("crossing") {
  using namespace dst;
  /*
    <---------+0+------->
    |                   |
    |                   |
    v                   v
    1+-------->11+----->2
    +                  ++
    |                  ||
    v                  +v
    21         22<-----+23
  */

  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(0,2), 
                                         std::make_pair(1,11), 
                                         std::make_pair(11,2),
                                         std::make_pair(1,21), 
                                         std::make_pair(2,22), 
                                         std::make_pair(2,23)
                                         };
  std::vector<double> weights {1,1, 0.2,0.2, 3,3,1};
  std::vector<int> terms {21,22,23};
  DST dt = DST(edges, weights, 0, terms);

  CHECK(dt.naive_alg().cost == 4+4+2-1);
  CHECK(dt.naive_alg().cost_sc == 4+4+2);
  auto tree = dt.level2_alg();
  CHECK(tree.cost_sc == 2+(1+3+0.2*2+3));
  CHECK(std::abs(tree.cost - tree.cost_sc) < EPSILON);
  CHECK(std::abs(tree.cost_trimmed() - (tree.cost_sc-0.2*2)) < EPSILON);
  CHECK((tree.terms_cov == std::unordered_set<int> {24,25,26}));
}


TEST_CASE("cycles") {
  using namespace dst;

  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(1,2), 
                                         std::make_pair(2,3), 
                                         std::make_pair(3,4), 
                                         std::make_pair(4,5),
                                         std::make_pair(5,0)};
  std::vector<double> weights {1,1,1,1,1,1};
  std::vector<int> terms {2,4};
  DST dt = DST(edges, weights, 0, terms);

  CHECK(dt.naive_alg().cost_sc == 2+4);
  CHECK(dt.naive_alg().cost == 2+2);
  auto tree = dt.level2_alg();
  CHECK(tree.cost == 2+2);
  // CHECK(tree.cost_sc == 2+4); // pick 2 and then 4
  CHECK(tree.cost_sc == 4); // pick 2 only
  CHECK((tree.terms_cov == std::unordered_set<int> {6,7}));
}


TEST_CASE("level3_alg") {
  using namespace dst;
  /*
   +------0-------+    
   |      +       |    
   1    +-3-+     2    
   |\   |   |    /|    
   | \ 11   12  | |    
   |  \/|    |\ / |    
   |  /\|    | \  |    
   | -  -    || \ |    
   21   22   23  -24  
  */

  // r - u - v - {t}
  // 4 + 2*1 + 4*1 = 10
  // (2 + 2*2) * 2 = 12

  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(0,2), 
                                         std::make_pair(0,3), 
                                         std::make_pair(1,21), 
                                         std::make_pair(1,22), 
                                         std::make_pair(2,23), 
                                         std::make_pair(2,24), 
                                         std::make_pair(3,11), 
                                         std::make_pair(3,12), 
                                         std::make_pair(11,21), 
                                         std::make_pair(11,22), 
                                         std::make_pair(12,23), 
                                         std::make_pair(12,24)};
  std::vector<double> weights {2,2,4, 2,2,2,2, 1.5,1.5, 1,1,1,1};
  std::vector<int> terms {21,22,23,24};
  DST dt = DST(edges, weights, 0, terms);

  auto &&tree_naive = dt.naive_alg();
  CHECK(tree_naive.cost_sc == (2+2) * 4);
  CHECK(tree_naive.cost == (2 + 2*2) * 2);

  auto &&tree2 = dt.level2_alg();
  CHECK(tree2.cost == (2 + 2*2) * 2);
  CHECK(tree2.cost_sc == (2 + 2*2) * 2);
  CHECK((tree2.terms_cov == std::unordered_set<int> {25,26,27,28}));

  auto &&tree3 = dt.level3_alg(0.99);
  CHECK(tree3.cost == 4 + 2*1.5 + 4*1);
  CHECK(tree3.cost_sc == 4 + 2*1.5 + 4*1);
  CHECK((tree3.terms_cov == std::unordered_set<int> {25,26,27,28}));
}


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
  std::vector<double> weights {1,1,1,2.5, 2,2,2, 1,1,1};
  std::vector<int> terms {11,12,13};
  DST dt = DST(edges, weights, 0, terms);

  CHECK(dt.naive_alg().cost == 1*3 + 2*3);
  CHECK(dt.naive_alg().cost_sc == 1*3 + 2*3);

  auto tree = dt.level2_alg();
  CHECK(tree.cost == 2.5 + 1*3);
  CHECK(tree.cost_sc == 2.5 + 1*3);
  CHECK((tree.terms_cov == std::unordered_set<int> {14,15,16}));

  tree = dt.level2_co_alg();
  CHECK(tree.cost == 2.5 + 1*3);
  CHECK(tree.cost_sc == 2.5 + 1*3);
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
  
  CHECK(dt.naive_alg().cost == 1+1+2-1);
  CHECK(dt.naive_alg().cost_sc == 1+1+2);
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

  auto p1 = dijkstra(dt.adj, dt.w, 0);
  auto dists = p1.first;
  auto trace = p1.second;
  CHECK(dists.at(0) == 0);
  CHECK(dists.at(1) == 1);
  CHECK(dists.at(2) == 1);
  CHECK(dists.at(3) == 1);
  CHECK(dists.at(4) == 2);
  CHECK(trace.at(4) == 3);

  auto p2 = dijkstra(dt.adj, dt.w, 4);
  dists = p2.first;
  trace = p2.second;
  CHECK (dists.find(0) == dists.end());
  CHECK (dists.find(1) == dists.end());
  CHECK (dists.find(2) == dists.end());
  CHECK (dists.at(3) == 1);
  CHECK (dists.at(4) == 0);

  auto p3 = dijkstra(dt.adj_r, dt.w, 4, true);
  dists = p3.first;
  trace = p3.second;
  CHECK (dists.at(0) == 2);
  CHECK (dists.find(1) == dists.end());
  CHECK (dists.find(2) == dists.end());
  CHECK (dists.at(3) == 1);
  CHECK (dists.at(4) == 0);
  CHECK (trace.at(4) == NONVERTEX);
  CHECK (trace.at(3) == 4);
  CHECK (trace.find(2) == trace.end());

  auto p4 = dijkstra(dt.adj, dt.w, 0, false, {3});
  CHECK (not has_key(p4.first, 4));
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

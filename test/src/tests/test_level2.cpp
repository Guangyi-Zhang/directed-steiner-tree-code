#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <tuple>

#include <fmt/ranges.h>
#include "dst/dst.hpp"
#include "dst/dijkstra.hpp"
#include "dst/utils.hpp"


TEST_CASE("fast_level2_alg") {
  using namespace dst;
  /*
    <-----+0+----->
    |             |
    |             |
    v             v
    1             2+------>
    +             +       |
    |             |       |
    |             |       v
    v----->11<----v       12
  */

  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(0,2), 
                                         std::make_pair(1,11), 
                                         std::make_pair(2,11), 
                                         std::make_pair(2,12), 
                                         };
  std::vector<double> weights {1,1, 99,1,1};
  int root = 0;
  std::vector<int> terms {11,12};
  DST dt = DST(edges, weights, root, terms);

  CHECK(dt.level1_alg()->cost == 3);
  CHECK(dt.level1_alg()->cost_sc == 4);
  CHECK(dt.level1_alg()->debuginfo.at("sssp_nodes_visited") == "7"); // 5+2 dummies

  auto tree = dt.level2_alg();

  CHECK(dt.level2_alg()->to_tree()->cost == 3);
  CHECK(dt.level2_alg()->cost_sc == 3);
  CHECK(dt.level2_alg()->debuginfo.at("sssp_nodes_visited") == "16"); // 7+(5+4)

  CHECK(dt.fast_level2_alg()->to_tree()->cost == 3);
  CHECK(dt.fast_level2_alg()->cost_sc == 3);
  CHECK(dt.fast_level2_alg()->debuginfo.at("sssp_nodes_visited") == "13"); // 7+2*(3)

  dt.fast_level2_alg()->to_tree()->print();
}


TEST_CASE("PartialTree3") {
  // test zero_drv()
  using namespace dst;

  int root {0}, v {1};
  double d_rv {1};
  int k = 5;
  PartialTree tree {root, v, d_rv};

  tree.add_term(11, 1); // density = 2/1
  CHECK(not tree.is_ready());
  tree.zero_drv();
  CHECK(not tree.is_ready()); // need at least 2 terminals
  CHECK(eq(tree.density(), 1));

  PartialTree tree2 {root, v, d_rv};
  tree2.add_term(11, 1); // density = 2/1
  tree2.add_term(12, 2); // density = 4/2
  tree2.add_term(13, 4); // density = 8/3
  CHECK(tree2.is_ready());
  CHECK(tree2.terms_after_ready.size() == 1);
  tree2.zero_drv();
  CHECK(tree2.is_ready());
  CHECK(eq(tree.density(), 1));
  CHECK(tree2.terms_after_ready.size() == 2);
  CHECK(tree2.terms_after_ready.front() == 12);
}


TEST_CASE("PartialTree2") {
  // test erase_and_reset()
  using namespace dst;

  int root {0}, v {1};
  double d_rv {1};
  int k = 5;
  PartialTree tree {root, v, d_rv};
  CHECK(not tree.is_ready());
  CHECK(tree.density() > 1e8);
  CHECK(tree.density_LB(k) < 0);

  tree.add_term(11, 1); // density = 2/1
  CHECK(eq(tree.density_LB(k), (1.0+k)/k));
  tree.add_term(12, 2); // density = 4/2
  CHECK(not tree.is_ready());
  CHECK(eq(tree.density_LB(k), (1.0+1+(k-1)*2)/k));
  tree.add_term(13, 4); // density = 8/3
  CHECK(tree.is_ready());
  CHECK(eq(tree.density(), 2));
  CHECK(eq(tree.density_LB(k), tree.density()));
  CHECK(tree.terms.size() == 2);
  CHECK(tree.terms_after_ready.size() == 1);

  tree.erase_and_reset({11}); // density = (1+2)/1
  CHECK(tree.is_ready());
  CHECK(eq(tree.density(), 3));
  CHECK(eq(tree.density_LB(k), tree.density()));
  CHECK(tree.terms.size() == 1);
  CHECK(tree.terms_after_ready.size() == 1); 
  CHECK(tree.terms_after_ready == std::list<int> {13}); 
}


TEST_CASE("PartialTree") {
  using namespace dst;

  int root {0}, v {1};
  double d_rv {1};
  int k = 5;
  PartialTree tree {root, v, d_rv};
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

  PartialTree tree2 {0.5};
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

  CHECK(dt.level1_alg()->cost == 4+4+2-1);
  CHECK(dt.level1_alg()->cost_sc == 4+4+2);

  // first pick 2->23, and then 1->{21,22}
  // no, now adaptively pick 2->23, then 2->22, then 1->21
  auto tree = dt.level2_alg();
  CHECK((tree->terms == std::unordered_set<int> {24,25,26}));
  CHECK(tree->cost_sc == 2+3+4);
  //CHECK(std::abs(tree->cost_sc - tree->to_tree()->cost - 0.2) < EPSILON); // include 1->11 but 11->2
  //CHECK(std::abs(tree->to_tree()->cost_trimmed() - (tree->cost_sc-0.2*2)) < EPSILON);
  CHECK(std::abs(tree->cost_sc - tree->to_tree()->cost) < EPSILON); // include 1->11 but 11->2
  CHECK(std::abs(tree->to_tree()->cost_trimmed() - (tree->cost_sc)) < EPSILON);
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

  CHECK(dt.level1_alg()->cost_sc == 2+4);
  CHECK(dt.level1_alg()->cost == 2+2);
  auto tree = dt.level2_alg();
  CHECK(tree->to_tree()->cost == 2+2);
  // CHECK(tree->cost_sc == 2+4); // pick 2 and then 4
  CHECK(tree->cost_sc == 4); // pick 2 only
  CHECK((tree->terms == std::unordered_set<int> {6,7}));
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

  CHECK(dt.level1_alg()->cost == 1*3 + 2*3);
  CHECK(dt.level1_alg()->cost_sc == 1*3 + 2*3);

  auto tree = dt.level2_alg();
  CHECK(tree->to_tree()->cost == 2.5 + 1*3);
  CHECK(tree->cost_sc == 2.5 + 1*3);
  CHECK((tree->terms == std::unordered_set<int> {14,15,16}));

  tree = dt.fast_level2_alg();
  CHECK(tree->to_tree()->cost == 2.5 + 1*3);
  CHECK(tree->cost_sc == 2.5 + 1*3);
  CHECK((tree->terms == std::unordered_set<int> {14,15,16}));
}


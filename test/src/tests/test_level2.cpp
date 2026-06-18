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


TEST_CASE("sssp_nodes_visited") {
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
  CHECK(dt.level1_alg()->debuginfo.at("sssp_nodes_visited") == "7"); // 5+2=7 dummies

  auto tree = dt.level2_alg();

  CHECK(dt.level2_alg()->to_tree()->cost == 3);
  CHECK(dt.level2_alg()->cost_sc == 3);
  CHECK(dt.level2_alg()->debuginfo.at("sssp_nodes_visited") == "16"); // 7+(5 for t11 + 4 for t12)

  auto fast2 = dt.fast_level2_alg();
  CHECK(fast2->to_tree()->cost == 3);
  CHECK(fast2->cost_sc == 3);
  CHECK(fast2->debuginfo.at("sssp_nodes_visited") == "14"); // 7+4+3

  fast2->to_tree()->print();
}


TEST_CASE("FenwickTree prefix sums and point updates") {
  using namespace dst;

  FenwickTree fw {6};
  CHECK(fw.size() == 6);
  CHECK(eq(fw.total_sum(), 0));
  CHECK(eq(fw.prefix_sum(0), 0));

  fw.add(0, 2);
  fw.add(2, 3.5);
  fw.add(5, 4);

  CHECK(eq(fw.at(0), 2));
  CHECK(eq(fw.at(2), 3.5));
  CHECK(eq(fw.prefix_sum(0), 2));
  CHECK(eq(fw.prefix_sum(1), 2));
  CHECK(eq(fw.prefix_sum(2), 5.5));
  CHECK(eq(fw.prefix_sum(5), 9.5));
  CHECK(eq(fw.total_sum(), 9.5));

  fw.set(2, 1.5);
  fw.add(4, 0.25);

  CHECK(eq(fw.at(2), 1.5));
  CHECK(eq(fw.prefix_sum(2), 3.5));
  CHECK(eq(fw.prefix_sum(4), 3.75));
  CHECK(eq(fw.total_sum(), 7.75));
}


TEST_CASE("FenwickTree select by cumulative count") {
  using namespace dst;

  FenwickTree counts {7};
  counts.add(1, 1);
  counts.add(3, 1);
  counts.add(4, 1);
  counts.add(6, 1);

  CHECK(eq(counts.total_sum(), 4));
  CHECK(counts.select_by_cumulative(0) == NONVERTEX);
  CHECK(counts.select_by_cumulative(1) == 1);
  CHECK(counts.select_by_cumulative(2) == 3);
  CHECK(counts.select_by_cumulative(3) == 4);
  CHECK(counts.select_by_cumulative(4) == 6);
  CHECK(counts.select_by_cumulative(5) == NONVERTEX);

  counts.add(3, -1);
  counts.set(6, 0);

  CHECK(eq(counts.total_sum(), 2));
  CHECK(counts.select_by_cumulative(1) == 1);
  CHECK(counts.select_by_cumulative(2) == 4);
  CHECK(counts.select_by_cumulative(3) == NONVERTEX);
}


TEST_CASE("AdaptiveLevel2CenterIndex builds threshold prefix and erases") {
  using namespace dst;

  AdaptiveLevel2CenterIndex idx {7, 4};
  idx.insert(101, 1);
  idx.insert(102, 2);
  idx.insert(103, 4);
  idx.insert(102, 2); // duplicate discoveries are ignored

  CHECK(idx.alive_discovered_count() == 3);
  CHECK(idx.best_prefix_length(1) == 2);

  auto bound = idx.discovered_bound(1);
  CHECK(bound.nterm == 2);
  CHECK(eq(bound.density, 2));

  auto tree = idx.partree(0, 1);
  CHECK(tree->v == 7);
  CHECK(eq(tree->cost_sc, 4));
  CHECK((tree->terms == std::unordered_set<int> {101, 102}));
  CHECK(eq(tree->distances_t.at(101), 1));
  CHECK(eq(tree->distances_t.at(102), 2));

  idx.erase({101, 999});
  CHECK(idx.alive_discovered_count() == 2);
  CHECK(idx.best_prefix_length(1) == 1);

  tree = idx.partree(0, 1);
  CHECK(eq(tree->cost_sc, 3));
  CHECK((tree->terms == std::unordered_set<int> {102}));

  CHECK(idx.best_prefix_length(10) == 2);
  tree = idx.partree(0, 10);
  CHECK(eq(tree->cost_sc, 16));
  CHECK((tree->terms == std::unordered_set<int> {102, 103}));
}


TEST_CASE("AdaptiveLevel2CenterIndex optimistic bound handles undiscovered tail") {
  using namespace dst;

  AdaptiveLevel2CenterIndex idx {5, 3};
  idx.insert(201, 2);

  auto discovered = idx.discovered_bound(10);
  CHECK(discovered.nterm == 1);
  CHECK(eq(discovered.density, 12));

  auto optimistic = idx.optimistic_bound(10, 3, 3, false);
  CHECK(optimistic.nterm == 3);
  CHECK(eq(optimistic.density, 6));

  auto exhausted = idx.optimistic_bound(10, 3, 3, true);
  CHECK(exhausted.nterm == 1);
  CHECK(eq(exhausted.density, 12));

  AdaptiveLevel2CenterIndex empty {6, 2};
  optimistic = empty.optimistic_bound(5, 2, 3, false);
  CHECK(optimistic.nterm == 2);
  CHECK(eq(optimistic.density, 5.5));
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
  CHECK(eq(tree.density_LB(k), d_rv/k));

  tree.add_term(11, 1); // density = 2/1
  CHECK(eq(tree.density_LB(k,1), (1.0+k)/k));
  tree.add_term(12, 2); // density = 4/2
  CHECK(not tree.is_ready());
  CHECK(eq(tree.density_LB(k,2), (1.0+1+(k-1)*2)/k));
  tree.add_term(13, 4); // density = 8/3
  CHECK(tree.is_ready());
  CHECK(eq(tree.density(), 2));
  CHECK(eq(tree.density_LB(k,4), tree.density()));
  CHECK(tree.terms.size() == 2);
  CHECK(tree.terms_after_ready.size() == 1);

  tree.erase_and_reset({11}); // density = (1+2)/1
  CHECK(tree.is_ready());
  CHECK(eq(tree.density(), 3));
  CHECK(eq(tree.density_LB(k,4), tree.density()));
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
  CHECK(eq(tree.density_LB(k), d_rv/k));

  tree.add_term(11, 1); // density = 2/1
  CHECK(eq(tree.density_LB(k,1), (1.0+k)/k));
  tree.add_term(12, 2); // density = 4/2
  CHECK(not tree.is_ready());
  CHECK(eq(tree.density_LB(k,2), (1.0+1+(k-1)*2)/k));
  tree.add_term(13, 3); // density = 7/3
  CHECK(tree.is_ready());
  CHECK(eq(tree.density(), 2));
  CHECK(eq(tree.density_LB(k,3), tree.density()));
  CHECK(tree.terms.size() == 2);

  tree.erase_and_reset({11, 12});
  CHECK(not tree.is_ready());
  CHECK(eq(tree.density(), 1+3));
  CHECK(eq(tree.density_LB(k,3), (1.0+k*3)/k));
  CHECK(tree.terms.size() == 1);

  tree.add_term(14, 4); // density = (1+3+4)/2
  CHECK(eq(tree.density_LB(k,4), (1.0+3+(k-1)*4)/k));
  tree.add_term(14, 5); // density = (1+3+4+5)/3
  CHECK(tree.is_ready());
  CHECK(eq(tree.density(), 4));
  CHECK(eq(tree.density_LB(k,4), tree.density()));
  CHECK(tree.terms.size() == 2);
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


TEST_CASE("fast_level2") {
  using namespace dst;

  // use an arc to prune all v's along a long path
  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(1,2), 
                                         std::make_pair(2,3), 
                                         std::make_pair(3,4), 
                                         std::make_pair(4,5), 
                                         
                                         std::make_pair(0,5), 
                                         };
  std::vector<double> weights {1,1,1,1,1,1};
  int root = 0;
  std::vector<int> terms {5};
  DST dt = DST(edges, weights, root, terms);
  DEBUG_fast_level2 = true;

  auto fast2 = dt.fast_level2_alg();
  //CHECK(dt.strtest == "5:1.000000,4:5.000000,3:5.000000,2:5.000000,1:5.000000,0:1.000000,greedy_by_pruning"); // v:denlb
  CHECK(dt.strtest == "5:1.000000,4:4.000000,3:3.000000,2:2.000000,1:1.000000,0:1.000000,greedy_by_pruning"); // prune when d_vt=0
  
  CHECK(eq(fast2->to_tree()->cost, 1));
  CHECK(eq(fast2->cost_sc, 1));
}


TEST_CASE("adaptive_level2_reuses_root_side_distance") {
  using namespace dst;

  /*
    <---+0+--->
          |(4)
          v
          1+--(1)-->11
          |(5)
          v
          2+--(1)-->21
          +
          |(1)
          v
         22

    After covering terminal 11 through vertex 1, vertex 2 is reachable from the
    current tree with marginal cost 5. fast_level2 keeps the original root-side
    distance to 2, while adaptive_level2 and full_adaptive_level2 pay only the
    marginal root-side distance for the second two-terminal bunch.
  */
  std::vector<std::pair<int,int>> edges {std::make_pair(0,1),
                                         std::make_pair(1,11),
                                         std::make_pair(1,2),
                                         std::make_pair(2,21),
                                         std::make_pair(2,22)};
  std::vector<double> weights {4,1,5,1,1};
  std::vector<int> terms {11,21,22};
  DST dt = DST(edges, weights, 0, terms);

  auto fast2 = dt.fast_level2_alg();
  auto adaptive2 = dt.adaptive_level2_alg();
  auto adaptive_fast2 = dt.adaptive_fast_level2_alg();
  auto full_adaptive2 = dt.full_adaptive_level2_alg();
  auto full_adaptive_fast2 = dt.full_adaptive_fast_level2_alg();

  CHECK(eq(fast2->cost_sc, 16));
  CHECK(eq(adaptive2->cost_sc, 12));
  CHECK(eq(adaptive_fast2->cost_sc, 12));
  CHECK(eq(full_adaptive2->cost_sc, 12));
  CHECK(eq(full_adaptive_fast2->cost_sc, 12));

  CHECK(eq(fast2->to_tree()->cost, 12));
  CHECK(eq(adaptive2->to_tree()->cost, 12));
  CHECK(eq(adaptive_fast2->to_tree()->cost, 12));
  CHECK(eq(full_adaptive2->to_tree()->cost, 12));
  CHECK(eq(full_adaptive_fast2->to_tree()->cost, 12));
}


TEST_CASE("full_adaptive_level2_reuses_terminal_side_residual_paths") {
  using namespace dst;

  /*
    <---+0+--(2)-->4
    |(8)  |
    v     +--(1)-->6+--(2)-->5
    2     |
          +--(4)-->5

    All three algorithms find the same final tree (cost 5): 0->6, 0->4, 6->5.
    fast_level2 and adaptive_level2 take three greedy bunches: {6}, then {5}
    through v=6, then {4} through v=0. The terminal-side distance 6->5 is
    always 2, so after {6} is covered, {5} alone through v=6 has density 2.
    Bunching {4,5} needs a shared middle vertex; through v=0 that costs
    d_4(0)+d_5(0)=2+3=5 (density 2.5). 
    full_adaptive_level2 takes two bunches: {6}, then {4,5} through v=0.
    Residual weights zero out 0->6, so the terminal-side distance d_5(0)
    drops to 2 via 5->6->0 and bunching {4,5} reaches density 2 as well.
  */
  std::vector<std::pair<int,int>> edges {std::make_pair(0,2),
                                         std::make_pair(0,4),
                                         std::make_pair(0,5),
                                         std::make_pair(0,6),
                                         std::make_pair(2,6),
                                         std::make_pair(3,1),
                                         std::make_pair(3,4),
                                         std::make_pair(5,0),
                                         std::make_pair(6,0),
                                         std::make_pair(6,3),
                                         std::make_pair(6,4),
                                         std::make_pair(6,5)};
  std::vector<double> weights {8,2,4,1,6,4,8,6,4,1,6,2};
  std::vector<int> terms {4,5,6};
  DST dt = DST(edges, weights, 0, terms);

  auto fast2 = dt.fast_level2_alg();
  auto adaptive2 = dt.adaptive_level2_alg();
  auto adaptive_fast2 = dt.adaptive_fast_level2_alg();
  auto full_adaptive2 = dt.full_adaptive_level2_alg();
  auto full_adaptive_fast2 = dt.full_adaptive_fast_level2_alg();

  CHECK(eq(fast2->cost_sc, 5));
  CHECK(eq(adaptive2->cost_sc, 5));
  CHECK(eq(adaptive_fast2->cost_sc, 5));
  CHECK(eq(full_adaptive2->cost_sc, 5));
  CHECK(eq(full_adaptive_fast2->cost_sc, 5));

  CHECK(fast2->subtrees.size() == 3);
  CHECK(adaptive2->subtrees.size() == 3);
  CHECK(adaptive_fast2->subtrees.size() == 3);
  CHECK(full_adaptive2->subtrees.size() == 2);
  CHECK(full_adaptive_fast2->subtrees.size() == 2);
  CHECK(full_adaptive2->subtrees.back().first->terms.size() == 2);
  CHECK(full_adaptive_fast2->subtrees.back().first->terms.size() == 2);
}


TEST_CASE("level2_alg") {
  using namespace dst;

  // large d(0,4) to fail level1
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

  tree = dt.adaptive_level2_alg();
  CHECK(tree->to_tree()->cost == 2.5 + 1*3);
  CHECK(tree->cost_sc == 2.5 + 1*3);
  CHECK((tree->terms == std::unordered_set<int> {14,15,16}));

  tree = dt.adaptive_fast_level2_alg();
  CHECK(tree->to_tree()->cost == 2.5 + 1*3);
  CHECK(tree->cost_sc == 2.5 + 1*3);
  CHECK((tree->terms == std::unordered_set<int> {14,15,16}));

  tree = dt.full_adaptive_level2_alg();
  CHECK(tree->to_tree()->cost == 2.5 + 1*3);
  CHECK(tree->cost_sc == 2.5 + 1*3);
  CHECK((tree->terms == std::unordered_set<int> {14,15,16}));

  tree = dt.full_adaptive_fast_level2_alg();
  CHECK(tree->to_tree()->cost == 2.5 + 1*3);
  CHECK(tree->cost_sc == 2.5 + 1*3);
  CHECK((tree->terms == std::unordered_set<int> {14,15,16}));
}

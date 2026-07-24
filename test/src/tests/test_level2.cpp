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
  std::vector<float> weights {1,1, 99,1,1};
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


TEST_CASE("level2 scored keys form a strict total order") {
  using namespace dst;

  struct Key { double density; int nterm; int center; };
  std::vector<Key> keys {
    {1.0, 2, 3},
    {1.0 + EPSILON * 0.6, 1, 2},
    {1.0 + EPSILON * 1.2, 3, 1},
    {1.0, 1, 4},
    {1.0, 1, 2}
  };
  auto less = [] (const Key &a, const Key &b) {
    return DST::level2_scored_key_less(
        a.density, a.nterm, a.center,
        b.density, b.nterm, b.center
    );
  };

  for (const auto &a : keys) {
    CHECK_FALSE(less(a, a));
    for (const auto &b : keys) {
      if (less(a, b))
        CHECK_FALSE(less(b, a));
      for (const auto &c : keys) {
        if (less(a, b) and less(b, c))
          CHECK(less(a, c));
      }
    }
  }

  std::sort(keys.begin(), keys.end(), less);
  CHECK(keys[0].density == 1.0);
  CHECK(keys[0].nterm == 1);
  CHECK(keys[0].center == 2);
  CHECK(keys[1].density == 1.0);
  CHECK(keys[1].nterm == 1);
  CHECK(keys[1].center == 4);
  CHECK(keys[2].density == 1.0);
  CHECK(keys[2].nterm == 2);
  CHECK(keys[3].density == 1.0 + EPSILON * 0.6);
  CHECK(keys[4].density == 1.0 + EPSILON * 1.2);
}


TEST_CASE("observed-prefix pool uses bounded ordered replacement") {
  using namespace dst;

  std::vector<DST::Level2ObservedPrefix> top;
  top.reserve(2);
  DST::Level2ObservedPoolStats stats;
  auto prefix = [] (double density, int nterm, int center, bool ready = false) {
    return DST::make_level2_observed_prefix(
        density, density * nterm, nterm, center, ready
    );
  };

  CHECK_FALSE(DST::consider_observed_prefix(top, prefix(0, 0, 9), 2, &stats));
  CHECK(DST::consider_observed_prefix(top, prefix(4, 2, 4), 2, &stats));
  CHECK(DST::consider_observed_prefix(top, prefix(3, 1, 3), 2, &stats));
  CHECK_FALSE(DST::consider_observed_prefix(top, prefix(5, 1, 5), 2, &stats));
  CHECK(top.size() == 2);
  CHECK(top[0].center == 3);
  CHECK(top[1].center == 4);

  CHECK(DST::consider_observed_prefix(top, prefix(2, 3, 4), 2, &stats));
  CHECK(top.size() == 2);
  CHECK(top[0].center == 4);
  CHECK(top[1].center == 3);

  CHECK_FALSE(DST::consider_observed_prefix(top, prefix(2, 3, 4, true), 2, &stats));
  CHECK(top[0].ready);
  CHECK_FALSE(DST::consider_observed_prefix(top, prefix(2.5, 3, 4), 2, &stats));
  CHECK(top[0].scored_density == 2);
  CHECK(stats.invariant_violations == 1);

  std::vector<DST::Level2ObservedPrefix> one;
  CHECK(DST::consider_observed_prefix(one, prefix(3, 2, 2), 1));
  CHECK(DST::consider_observed_prefix(one, prefix(2, 2, 1), 1));
  CHECK(one.size() == 1);
  CHECK(one.front().center == 1);
}


TEST_CASE("bounded observed-prefix pool matches a sorted reference stream") {
  using namespace dst;

  constexpr int rerank_L = 5;
  std::vector<DST::Level2ObservedPrefix> bounded;
  std::vector<DST::Level2ObservedPrefix> reference;
  bounded.reserve(rerank_L);

  for (int i = 0; i < 200; i++) {
    int center = (i * 7) % 23;
    int observation = i / 23;
    double density = 20.0 + center * 0.1 - observation * 0.01;
    int nterm = 1 + observation;
    auto cand = DST::make_level2_observed_prefix(
        density, density * nterm, nterm, center, (i % 3) == 0
    );
    DST::consider_observed_prefix(bounded, cand, rerank_L);

    auto old = std::find_if(reference.begin(), reference.end(), [&] (const auto &item) {
      return item.center == center;
    });
    if (old != reference.end()) {
      if (DST::level2_observed_prefix_less(*old, cand)) {
        old->ready = cand.ready;
      } else {
        reference.erase(old);
        reference.push_back(cand);
      }
    } else {
      reference.push_back(cand);
    }
    std::sort(
        reference.begin(),
        reference.end(),
        DST::level2_observed_prefix_less
    );
    if (static_cast<int>(reference.size()) > rerank_L)
      reference.resize(rerank_L);

    REQUIRE(bounded.size() == reference.size());
    for (size_t j = 0; j < bounded.size(); j++) {
      CHECK(bounded[j].center == reference[j].center);
      CHECK(bounded[j].nterm == reference[j].nterm);
      CHECK(bounded[j].scored_density == reference[j].scored_density);
      CHECK(bounded[j].scored_cost == reference[j].scored_cost);
    }
  }
}


TEST_CASE("level2 center indexes reconstruct observed prefixes") {
  using namespace dst;

  Level2CenterIndex idx(0, 7, 10.0);
  CHECK(idx.add_term(101, 1));
  double snapshot_cost = idx.tree.cost_sc;
  CHECK(idx.add_term(102, 2));
  auto stale = idx.partree_prefix(1, false);
  CHECK(stale->v == 7);
  CHECK(stale->terms == std::unordered_set<int> {101});
  CHECK(eq(stale->cost_sc, snapshot_cost));

  AdaptiveLevel2CenterIndex adaptive {8, 5};
  adaptive.insert(201, 1);
  adaptive.insert(202, 1);
  adaptive.insert(203, 2);
  adaptive.discovered_bound(10);
  adaptive.erase({201});
  auto ranked = adaptive.partree(0, 10, 2);
  CHECK(ranked->terms == std::unordered_set<int> {202, 203});
  CHECK(ranked->terms_before_ready.front() == 202);
  CHECK(ranked->terms_before_ready.back() == 203);
}


TEST_CASE("level2 center indexes report selected-prefix changes") {
  using namespace dst;

  Level2CenterIndex idx(0, 1, 1.0);
  CHECK(idx.add_term(11, 1));
  CHECK(idx.add_term(12, 3)); // readiness changes even though the prefix does not
  CHECK_FALSE(idx.add_term(13, 5));
  CHECK_FALSE(idx.add_term(12, 3));

  AdaptiveLevel2CenterIndex adaptive {5, 4};
  adaptive.insert(21, 1);
  adaptive.insert(22, 3);
  adaptive.discovered_bound(10);
  CHECK(adaptive.insert(23, 5));
  CHECK_FALSE(adaptive.insert(24, 20));
  CHECK_FALSE(adaptive.insert(24, 20));
}


TEST_CASE("AdaptiveLevel2CenterIndex builds threshold prefix and erases") {
  using namespace dst;

  AdaptiveLevel2CenterIndex idx {7, 4};
  idx.insert(101, 1);
  idx.insert(102, 2);
  idx.insert(103, 4);
  idx.insert(102, 2); // duplicate discoveries are ignored

  CHECK(idx.alive_discovered_count() == 3);
  CHECK(idx.best_prefix_length(1) == 1);

  auto bound = idx.discovered_bound(1);
  CHECK(bound.nterm == 1);
  CHECK(eq(bound.density, 2));

  auto tree = idx.partree(0, 1);
  CHECK(tree->v == 7);
  CHECK(eq(tree->cost_sc, 2));
  CHECK((tree->terms == std::unordered_set<int> {101}));
  CHECK(eq(tree->distances_t.at(101), 1));

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


TEST_CASE("AdaptiveLevel2CenterIndex certificate handles undiscovered tail") {
  using namespace dst;

  AdaptiveLevel2CenterIndex idx {5, 3};
  idx.insert(201, 2);

  auto discovered = idx.discovered_bound(10);
  CHECK(discovered.nterm == 1);
  CHECK(eq(discovered.density, 12));

  // Full-tail case: unseen terminals are cheap enough to possibly improve the
  // prefix, so the certificate uses the full-tail lower bound.
  auto optimistic = idx.optimistic_bound(10, 3, 3, false);
  CHECK(optimistic.nterm == 3);
  CHECK(eq(optimistic.density, 6));

  // Exhausted terminal search: no unseen terminals remain, so the discovered
  // best prefix is exact.
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
  PartialTree tree {root, v, d_rv};

  tree.add_term(11, 1); // density = 2/1
  CHECK(not tree.is_ready());
  tree.zero_drv();
  CHECK(not tree.is_ready()); // need at least 2 terminals
  CHECK(eq(tree.density(), 1));

  PartialTree tree2 {root, v, d_rv};
  tree2.add_term(11, 1); // density = 2/1
  tree2.add_term(12, 2); // density = 4/2, equal to current density, so stop
  tree2.add_term(13, 4); // density = 8/3
  CHECK(tree2.is_ready());
  CHECK(tree2.terms_after_ready.size() == 2);
  tree2.zero_drv();
  CHECK(tree2.is_ready());
  CHECK(eq(tree2.density(), 1));
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
  tree.add_term(12, 2); // density = 4/2, equal to current density, so stop
  CHECK(tree.is_ready());
  CHECK(eq(tree.density_LB(k,2), (1.0+1+(k-1)*2)/k));
  tree.add_term(13, 4); // density = 8/3
  CHECK(tree.is_ready());
  CHECK(eq(tree.density(), 2));
  CHECK(eq(tree.density_LB(k,4), tree.density()));
  CHECK(tree.terms.size() == 1);
  CHECK(tree.terms_after_ready.size() == 2);

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
  tree.add_term(12, 2); // density = 4/2, equal to current density, so stop
  CHECK(tree.is_ready());
  CHECK(eq(tree.density_LB(k,2), (1.0+1+(k-1)*2)/k));
  tree.add_term(13, 3); // density = 7/3
  CHECK(tree.is_ready());
  CHECK(eq(tree.density(), 2));
  CHECK(eq(tree.density_LB(k,3), tree.density()));
  CHECK(tree.terms.size() == 1);

  tree.erase_and_reset({11, 12});
  CHECK(not tree.is_ready());
  CHECK(eq(tree.density(), 1+3));
  CHECK(eq(tree.density_LB(k,3), (1.0+k*3)/k));
  CHECK(tree.terms.size() == 1);

  tree.add_term(14, 4); // density = (1+3+4)/2
  CHECK(tree.is_ready());
  CHECK(eq(tree.density_LB(k,4), tree.density()));
  tree.add_term(14, 5); // density = (1+3+4+5)/3
  CHECK(tree.is_ready());
  CHECK(eq(tree.density(), 4));
  CHECK(eq(tree.density_LB(k,4), tree.density()));
  CHECK(tree.terms.size() == 1);
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
  std::vector<float> weights {1,1, 0.2,0.2, 3,3,1};
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
  std::vector<float> weights {1,1,1,1,1,1};
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
  std::vector<float> weights {1,1,1,1,1,1};
  int root = 0;
  std::vector<int> terms {5};
  DST dt = DST(edges, weights, root, terms);
  DEBUG_fast_level2 = true;

  auto fast2 = dt.fast_level2_alg();
  //CHECK(dt.strtest == "5:1.000000,4:5.000000,3:5.000000,2:5.000000,1:5.000000,0:1.000000,greedy_by_pruning"); // v:denlb
  CHECK(dt.strtest == "0:1.000000,1:2.000000,2:3.000000,3:4.000000,4:5.000000,5:1.000000,greedy_by_pruning"); // prune when d_vt=0
  
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
  std::vector<float> weights {4,1,5,1,1};
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
    Under the prefer-fewer-terminals tie rule, full_adaptive_level2 still takes
    three bunches: {6}, then one of {4,5}, then the other.
    Residual weights zero out 0->6, so the terminal-side distance d_5(0)
    drops to 2 via 5->6->0, but bunching {4,5} only ties density 2.
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
  std::vector<float> weights {8,2,4,1,6,4,8,6,4,1,6,2};
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
  CHECK(full_adaptive2->subtrees.size() == 3);
  CHECK(full_adaptive_fast2->subtrees.size() == 3);
  CHECK(full_adaptive2->subtrees.back().first->terms.size() == 1);
  CHECK(full_adaptive_fast2->subtrees.back().first->terms.size() == 1);
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
  std::vector<float> weights {1,1,1,2.5, 2,2,2, 1,1,1};
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


TEST_CASE("partree_realized_delta_scratch matches copy evaluator") {
  using namespace dst;

  /*
          0
          |(5)
          v
          1+--(3)-->2+--(5)-->11
          |
          +(4)-->3+--(6)-->12

    The base tree already paid 0->1.  The candidate rooted at 1 should only
    charge terminal-side arcs, and the scratch evaluator should match the
    copy-based evaluator exactly.
  */
  std::vector<std::pair<int,int>> edges {
    {0,1}, {1,2}, {2,11}, {1,3}, {3,12}
  };
  std::vector<float> weights {5,3,5,4,6};
  std::vector<int> terms {11,12};
  DST dt = DST(edges, weights, 0, terms);

  auto [dists_r, trace_r] = dijkstra(dt.graph.fwd_offset, dt.graph.fwd_target,
                                     dt.graph.fwd_weight, dt.root);
  auto trace_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,int>>>>();
  auto dists_t = std::make_shared<std::unordered_map<int, std::shared_ptr<std::unordered_map<int,float>>>>();
  for (auto t: dt.terms_dm) {
    auto [dists_, trace_] = dijkstra(dt.graph.rev_offset, dt.graph.rev_target,
                                     dt.graph.rev_weight, t);
    (*dists_t)[t] = dists_;
    (*trace_t)[t] = trace_;
  }

  auto base = std::make_shared<Tree> (dt.root);
  base->add_arc(std::make_pair(dt.root, 1), dists_r->at(1), trace_r);

  std::unordered_set<int> terms_left(dt.terms_dm.begin(), dt.terms_dm.end());
  auto cand = dt.level2_through_v(dt.root, 1, dists_r->at(1), dists_t, terms_left);

  double copy_delta = dt.partree_realized_delta(base, cand, trace_r, trace_t);
  double scratch_delta = dt.partree_realized_delta_scratch(base, cand, trace_r, trace_t);
  CHECK(eq(copy_delta, scratch_delta));
  CHECK(eq(scratch_delta, 18));

  double aborted = dt.partree_realized_delta_scratch(base, cand, trace_r, trace_t, 2);
  CHECK(aborted >= 2);
  CHECK(aborted <= scratch_delta);
}


TEST_CASE("adaptive level2 rerank_L one remains supported") {
  using namespace dst;

  std::vector<std::pair<int,int>> edges {
    {0,1}, {0,2}, {0,3}, {0,4},
    {1,11}, {2,12}, {3,13},
    {4,11}, {4,12}, {4,13}
  };
  std::vector<float> weights {1,1,1,2.5, 2,2,2, 1,1,1};
  std::vector<int> terms {11,12,13};
  DST dt = DST(edges, weights, 0, terms);

  auto adaptive2 = dt.adaptive_level2_alg(1);
  auto adaptive_fast2 = dt.adaptive_fast_level2_alg(1);
  auto full_adaptive2 = dt.full_adaptive_level2_alg(1);
  auto full_adaptive_fast2 = dt.full_adaptive_fast_level2_alg(1);

  CHECK(eq(adaptive2->to_tree()->cost, 5.5));
  CHECK(eq(adaptive_fast2->to_tree()->cost, 5.5));
  CHECK(eq(full_adaptive2->to_tree()->cost, 5.5));
  CHECK(eq(full_adaptive_fast2->to_tree()->cost, 5.5));

  CHECK(adaptive2->debuginfo.at("rerank_L") == "1");
  CHECK(adaptive_fast2->debuginfo.at("rerank_L") == "1");
  CHECK(full_adaptive2->debuginfo.at("rerank_L") == "1");
  CHECK(full_adaptive_fast2->debuginfo.at("rerank_L") == "1");
  CHECK(adaptive_fast2->debuginfo.at("rerank_mode") == "observed_prefix");
  CHECK(full_adaptive_fast2->debuginfo.at("rerank_mode") == "observed_prefix");
  CHECK(adaptive_fast2->debuginfo.at("certificate_L") == "1");
  CHECK(full_adaptive_fast2->debuginfo.at("certificate_L") == "1");
}


TEST_CASE("adaptive level2 APIs reject nonpositive rerank_L") {
  using namespace dst;

  std::vector<std::pair<int,int>> edges {{0, 1}, {1, 2}};
  std::vector<float> weights {1, 1};
  DST dt(edges, weights, 0, {2});

  CHECK_THROWS_AS(dt.adaptive_level2_alg(0), std::invalid_argument);
  CHECK_THROWS_AS(dt.adaptive_fast_level2_alg(-1), std::invalid_argument);
  CHECK_THROWS_AS(dt.full_adaptive_level2_alg(0), std::invalid_argument);
  CHECK_THROWS_AS(dt.full_adaptive_fast_level2_alg(-1), std::invalid_argument);
}


TEST_CASE("fast adaptive observed pool materializes only the final shortlist") {
  using namespace dst;

  std::vector<std::pair<int,int>> edges {
    {0,1}, {0,2}, {0,3}, {0,4},
    {1,11}, {2,12}, {3,13},
    {4,11}, {4,12}, {4,13}
  };
  std::vector<float> weights {1,1,1,2.5, 2,2,2, 1,1,1};
  DST dt(edges, weights, 0, {11,12,13});

  auto adaptive = dt.adaptive_fast_level2_alg(2);
  auto full = dt.full_adaptive_fast_level2_alg(2);
  for (const auto &result : {adaptive, full}) {
    auto considered = std::stoll(result->debuginfo.at("observed_pool_consider_calls"));
    auto materialized = std::stoll(result->debuginfo.at("observed_pool_materializations"));
    auto evaluated = std::stoll(result->debuginfo.at("rerank_candidates_evaluated"));
    CHECK(considered >= materialized);
    CHECK(materialized == evaluated);
    CHECK(materialized <= 2 * static_cast<long long>(result->subtrees.size()));
    CHECK(result->debuginfo.at("observed_pool_invariant_violations") == "0");
  }
}


TEST_CASE("level2_threshold returns q * next_weight - prefix_weight") {
  using namespace dst;
  CHECK(eq(level2_threshold(2, 5.0, 3.0), 7.0));
  CHECK(eq(level2_threshold(0, 1.0, 0.0), 0.0));
  CHECK(eq(level2_threshold(1, 4.0, 2.0), 2.0));
}


TEST_CASE("best_level2_prefix_length with plain vector thresholds") {
  using namespace dst;

  // Terminals with weights [1, 2, 4]: prefix_sums [1, 3, 7]
  // theta(1) = 1*2 - 1 = 1, theta(2) = 2*4 - 3 = 5
  std::vector<float> weights {1, 2, 4};
  std::vector<float> prefix_sums {1, 3, 7};
  auto thr = [&](int q) {
    return level2_threshold(q, weights[q], prefix_sums[q-1]);
  };

  // Equal thresholds keep the shorter prefix.
  CHECK(best_level2_prefix_length(3, 0, thr) == 1);
  CHECK(best_level2_prefix_length(3, 1, thr) == 1);
  CHECK(best_level2_prefix_length(3, 5, thr) == 2);
  CHECK(best_level2_prefix_length(3, 100, thr) == 3);

  auto dummy = [](int) { return 0.0; };
  CHECK(best_level2_prefix_length(0, 10, dummy) == 0);
  CHECK(best_level2_prefix_length(1, 10, dummy) == 1);
}


TEST_CASE("Level2CenterIndex ready center certifies the discovered prefix") {
  using namespace dst;
  Level2CenterIndex idx(0, 1, 1.0);
  idx.add_term(101, 1);
  idx.add_term(102, 2);
  idx.add_term(103, 4);

  // After adding 102, tree becomes ready with prefix {101}, cost_sc=2 because
  // the second terminal only ties the current density. The discovered rejected terminal proves the
  // certificate uses the current prefix density regardless of the tail.
  auto bound = idx.optimistic_bound(1.0, 3, 0, true);
  CHECK(bound.nterm == 1);
  CHECK(eq(bound.density, 2.0));
}


TEST_CASE("Level2CenterIndex certificate chooses current prefix when tail expensive") {
  using namespace dst;
  Level2CenterIndex idx(0, 1, 1.0);
  idx.add_term(11, 1);
  // Not ready, m=1, cost_sc=2. k=3, terminal_lb=10.
  // current prefix density = 2/1 = 2
  // full tail density = (2 + 2*10)/3 = 22/3 ~ 7.33
  // Expensive frontier: every unseen terminal (>= 10) is too expensive to
  // improve the current prefix, so the certificate uses the current prefix.
  auto bound = idx.optimistic_bound(1.0, 3, 10, false);
  CHECK(bound.nterm == 1);
  CHECK(eq(bound.density, 2.0));
}


TEST_CASE("Level2CenterIndex certificate chooses full tail when tail cheap") {
  using namespace dst;
  Level2CenterIndex idx(0, 1, 5.0);
  idx.add_term(11, 1);
  // Not ready, m=1, cost_sc=6. k=3, terminal_lb=1.
  // current prefix density = 6/1 = 6
  // full tail density = (6 + 2*1)/3 = 8/3 ~ 2.67
  // Cheap tail: unseen terminals could still improve the prefix, so the
  // certificate uses the full-tail lower bound (lem:level2-density-LB).
  auto bound = idx.optimistic_bound(5.0, 3, 1, false);
  CHECK(bound.nterm == 3);
  CHECK(eq(bound.density, 8.0 / 3));
}


TEST_CASE("Level2CenterIndex erase_and_reset preserves resume") {
  using namespace dst;
  Level2CenterIndex idx(0, 1, 1.0);
  idx.add_term(11, 1);
  idx.add_term(12, 2);
  idx.add_term(13, 4);
  CHECK(idx.tree.ready);

  idx.erase_and_reset({11});
  CHECK(idx.tree.terms.count(12) == 1); // contains 12

  // After erasing 11: prefix is {12}, cost_sc=3, ready. The discovered rejected
  // terminal 13 certifies the prefix exactly, so the certificate uses the
  // current prefix density regardless of the tail (terminals exhausted here).
  auto bound = idx.optimistic_bound(1.0, 4, 1.0, false);
  CHECK(bound.nterm == 1);
  CHECK(eq(bound.density, 3.0));
}


TEST_CASE("Level2CenterIndex zero_drv preserves certificate state") {
  using namespace dst;
  Level2CenterIndex idx(0, 1, 5.0);
  idx.add_term(11, 1);
  idx.add_term(12, 2);

  idx.zero_drv();

  // After zero_drv with d_rv=5: cost_sc resets to 0 and re-adds terminals.
  // Re-add 11 (d=1): cost_sc=1, density=1. Re-add 12 (d=2): new_den=1.5 > 1,
  // so tree becomes ready with prefix {11}, cost_sc=1.
  auto bound = idx.optimistic_bound(0.0, 2, 0, true);
  CHECK(bound.nterm == 1);
  CHECK(eq(bound.density, 1.0));

  // Same result regardless of d_rv parameter (d_rv is baked into cost_sc).
  auto bound2 = idx.optimistic_bound(1.0, 2, 0, true);
  CHECK(bound2.nterm == 1);
  CHECK(eq(bound2.density, 1.0));
}


TEST_CASE("Covered terminal in terms_after_ready is not re-added") {
  using namespace dst;
  // d_rv=1, terminals: 11(d=1), 12(d=3), 13(d=5)
  // add_term(11,1): cost_sc=2, terms={11}, density=2, not ready
  // add_term(12,3): new_den=(2+3)/2=2.5, lq(2,2.5)=true → ready, terms_after_ready=[12]
  // add_term(13,5): ready → terms_after_ready=[12,13]
  Level2CenterIndex idx(0, 1, 1.0);
  idx.add_term(11, 1);
  idx.add_term(12, 3);
  idx.add_term(13, 5);
  CHECK(idx.tree.ready);
  CHECK(idx.tree.terms.count(13) == 0);

  // Cover terminal 13. During replay, 12 is re-added and makes tree ready
  // before 13 is reached, so 13 stays in terms_after_ready.
  idx.erase_and_reset({13});
  CHECK(idx.tree.ready);
  CHECK(idx.tree.terms.count(13) == 0);

  // Now cover terminals 11 and 12. Replay reaches 13 but it must be skipped
  // because it was covered in the earlier batch.
  idx.erase_and_reset({11, 12});
  CHECK(idx.tree.terms.count(13) == 0);
  CHECK(idx.tree.terms.count(11) == 0);
  CHECK(idx.tree.terms.count(12) == 0);
}


TEST_CASE("Level2CenterIndex add_term ignores dead terminals") {
  using namespace dst;
  Level2CenterIndex idx(0, 1, 1.0);
  idx.add_term(11, 1);
  idx.add_term(12, 2);

  idx.erase_and_reset({12});

  // Re-adding a dead terminal must be a no-op.
  idx.add_term(12, 2);
  CHECK(idx.tree.terms.count(12) == 0);
}


TEST_CASE("AdaptiveLevel2CenterIndex certificate uses full tail when all discovered") {
  using namespace dst;
  AdaptiveLevel2CenterIndex idx {5, 4};
  idx.insert(201, 2);
  idx.insert(202, 3);

  // cached prefix at d_rv=10 uses both discovered terminals (nterm=m=2).
  // total_alive=4, unknown=2, terminal_lb=1.
  // current_density = (10 + 2 + 3) / 2 = 7.5
  // full_density = (10 + 2 + 3 + 2*1) / 4 = 4.25
  // Full tail is cheaper, so the certificate uses the full-tail lower bound.
  idx.discovered_bound(10);
  auto bound = idx.optimistic_bound(10, 4, 1, false);
  CHECK(bound.nterm == 4);
  CHECK(eq(bound.density, (10.0 + 2 + 3 + 2*1.0) / 4));
}


TEST_CASE("AdaptiveLevel2CenterIndex certificate is exact when search exhausted") {
  using namespace dst;
  AdaptiveLevel2CenterIndex idx {5, 3};
  idx.insert(201, 2);
  idx.insert(202, 3);
  idx.insert(203, 5);

  // Exhausted terminal search: no unseen terminals remain, so the discovered
  // best prefix is exact. At d_rv=10 every discovered terminal is taken:
  // theta(1)=1*3-2=1<=10, theta(2)=2*5-5=5<=10, so the prefix uses all three.
  idx.discovered_bound(10);
  auto bound = idx.optimistic_bound(10, 3, 0, true);
  CHECK(bound.nterm == 3);
  CHECK(eq(bound.density, (10.0 + 2 + 3 + 5) / 3));
}


TEST_CASE("fast_level2 and adaptive_fast_level2 single behavior") {
  using namespace dst;

  /*
              11   12   13
               ^    ^    ^
            (2)|    |    |(2)
               |    |    |
               1    2    3
               ^    ^    ^
            (1)|    |    |(1)
               |    |    |
              <---+0+---+
                    |
                (2.5)
                    v
                    4+--(1)-->11
                    +--(1)-->12
                    +--(1)-->13

    Optimal tree: 0->4 (2.5) plus 4->{11,12,13} (1 each), cost 5.5.
  */
  std::vector<std::pair<int,int>> edges {
    {0,1}, {0,2}, {0,3}, {0,4},
    {1,11}, {2,12}, {3,13},
    {4,11}, {4,12}, {4,13}
  };
  std::vector<float> weights {1,1,1,2.5, 2,2,2, 1,1,1};
  std::vector<int> terms {11,12,13};
  DST dt = DST(edges, weights, 0, terms);

  auto f2 = dt.fast_level2_alg();
  auto af2 = dt.adaptive_fast_level2_alg();

  CHECK(eq(f2->to_tree()->cost, 5.5));
  CHECK(eq(af2->to_tree()->cost, 5.5));
  CHECK(eq(f2->cost_sc, 5.5));
  CHECK(eq(af2->cost_sc, 5.5));
}


TEST_CASE("single behavior on adaptive root distance graph") {
  using namespace dst;

  /*
          0
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

    fast_level2 keeps root-side distance 9 to vertex 2, so cost_sc=16.
    adaptive_fast_level2 pays marginal cost 5 after covering 11, so cost_sc=12.
    Both build the same final tree (cost 12).
  */
  std::vector<std::pair<int,int>> edges {
    {0,1}, {1,11}, {1,2}, {2,21}, {2,22}
  };
  std::vector<float> weights {4,1,5,1,1};
  std::vector<int> terms {11,21,22};
  DST dt = DST(edges, weights, 0, terms);

  auto f2 = dt.fast_level2_alg();
  auto af2 = dt.adaptive_fast_level2_alg();

  CHECK(eq(f2->cost_sc, 16));
  CHECK(eq(af2->cost_sc, 12));
  CHECK(eq(f2->to_tree()->cost, 12));
  CHECK(eq(af2->to_tree()->cost, 12));
}


TEST_CASE("AdaptiveLevel2CenterIndex cached discovered bound matches binary search") {
  using namespace dst;

  AdaptiveLevel2CenterIndex idx {5, 5};
  idx.insert(201, 1);
  idx.insert(202, 3);
  idx.insert(203, 6);

  // theta(1) = 1*3 - 1 = 2 <= 4, theta(2) = 2*6 - 4 = 8 > 4
  // Best prefix length = 2, density = (4 + 1 + 3) / 2 = 4.
  auto bound = idx.discovered_bound(4.0);
  CHECK(bound.nterm == 2);
  CHECK(eq(bound.density, 4.0));

  // Second call with same d_rv uses cache.
  auto bound2 = idx.discovered_bound(4.0);
  CHECK(bound2.nterm == 2);
  CHECK(eq(bound2.density, 4.0));
}


TEST_CASE("AdaptiveLevel2CenterIndex cached prefix stops early ignores unknown") {
  using namespace dst;

  AdaptiveLevel2CenterIndex idx {5, 5};
  idx.insert(201, 1);
  idx.insert(202, 3);
  idx.insert(203, 6);

  // Refresh cache at d_rv = 4. Best prefix = 2 (stops at 203).
  idx.discovered_bound(4.0);

  // Unknown terminals at terminal_lb = 7 >= w_3 = 6:
  // prefix stops early at nterm=2, unknown tail cannot help.
  auto bound = idx.optimistic_bound(4.0, 5, 7.0, false);
  CHECK(bound.nterm == 2);
  CHECK(eq(bound.density, 4.0));
}


TEST_CASE("AdaptiveLevel2CenterIndex cached prefix full, cheap tail wins") {
  using namespace dst;

  AdaptiveLevel2CenterIndex idx {5, 3};
  idx.insert(201, 2);

  // Refresh cache at d_rv = 10. Best prefix = 1 = m.
  idx.discovered_bound(10.0);

  // total_alive=3, unknown=2, terminal_lb=3.
  // current_density = (10+2)/1 = 12
  // full_density = (10+2+2*3)/3 = 6
  // full tail wins.
  auto bound = idx.optimistic_bound(10.0, 3, 3.0, false);
  CHECK(bound.nterm == 3);
  CHECK(eq(bound.density, 6.0));
}


TEST_CASE("AdaptiveLevel2CenterIndex cached prefix full, expensive tail loses") {
  using namespace dst;

  AdaptiveLevel2CenterIndex idx {5, 3};
  idx.insert(201, 1);

  // Refresh cache at d_rv = 1. Best prefix = 1 = m.
  idx.discovered_bound(1.0);

  // total_alive=3, unknown=2, terminal_lb=10.
  // current_density = (1+1)/1 = 2
  // full_density = (1+1+2*10)/3 = 22/3 ~ 7.33
  // current prefix wins.
  auto bound = idx.optimistic_bound(1.0, 3, 10.0, false);
  CHECK(bound.nterm == 1);
  CHECK(eq(bound.density, 2.0));
}


TEST_CASE("AdaptiveLevel2CenterIndex insert extends cached prefix") {
  using namespace dst;

  AdaptiveLevel2CenterIndex idx {5, 5};
  idx.insert(201, 1);
  idx.insert(202, 3);

  // Refresh cache at d_rv = 10. Both terminals in prefix (nterm=2=m).
  auto bound = idx.discovered_bound(10.0);
  CHECK(bound.nterm == 2);
  CHECK(eq(bound.density, 7.0));

  // Insert terminal 203 with d=5. theta(2) = 2*5 - 4 = 6 <= 10, so it extends.
  idx.insert(203, 5);
  auto bound2 = idx.discovered_bound(10.0);
  CHECK(bound2.nterm == 3);
  CHECK(eq(bound2.density, (10.0 + 1 + 3 + 5) / 3));

  // Insert terminal 204 with d=20. theta(3) = 3*20 - 9 = 51 > 10, does not extend.
  idx.insert(204, 20);
  auto bound3 = idx.discovered_bound(10.0);
  CHECK(bound3.nterm == 3);
  CHECK(eq(bound3.density, (10.0 + 1 + 3 + 5) / 3));
}


TEST_CASE("AdaptiveLevel2CenterIndex erase invalidates cache, refresh restores") {
  using namespace dst;

  AdaptiveLevel2CenterIndex idx {5, 5};
  idx.insert(201, 1);
  idx.insert(202, 3);
  idx.insert(203, 6);

  idx.discovered_bound(4.0);

  // Erase terminal 201 (rank 0).
  idx.erase({201});

  // Alive terminals: 202 (d=3), 203 (d=6).
  // theta(1) = 1*6 - 3 = 3 <= 4, so nterm = 2.
  // density = (4 + 3 + 6) / 2 = 6.5
  auto bound = idx.discovered_bound(4.0);
  CHECK(bound.nterm == 2);
  CHECK(eq(bound.density, 6.5));
}


TEST_CASE("AdaptiveLevel2CenterIndex certificate matches expected on multiple configs") {
  using namespace dst;

  {
    // Config 1: cached prefix stops early, unknown tail ignored
    AdaptiveLevel2CenterIndex idx {5, 5};
    idx.insert(100, 1);
    idx.insert(101, 3);
    idx.insert(102, 6);
    idx.discovered_bound(4.0);
    auto r = idx.optimistic_bound(4.0, 5, 7.0, false);
    CHECK(r.nterm == 2);
    CHECK(eq(r.density, 4.0));
  }

  {
    // Config 2: all discovered, cheap tail wins
    AdaptiveLevel2CenterIndex idx {5, 3};
    idx.insert(100, 2);
    idx.discovered_bound(10.0);
    auto r = idx.optimistic_bound(10.0, 3, 3.0, false);
    CHECK(r.nterm == 3);
    CHECK(eq(r.density, 6.0));
  }

  {
    // Config 3: all discovered, expensive tail loses
    AdaptiveLevel2CenterIndex idx {5, 3};
    idx.insert(100, 1);
    idx.discovered_bound(1.0);
    auto r = idx.optimistic_bound(1.0, 3, 10.0, false);
    CHECK(r.nterm == 1);
    CHECK(eq(r.density, 2.0));
  }

  {
    // Config 4: exhausted search, no unknowns
    AdaptiveLevel2CenterIndex idx {5, 3};
    idx.insert(100, 1);
    idx.insert(101, 2);
    idx.insert(102, 4);
    idx.discovered_bound(1.0);
    auto r = idx.optimistic_bound(1.0, 3, 0, true);
    CHECK(r.nterm == 1);
    CHECK(eq(r.density, 2.0));
  }

  {
    // Config 5: no terminals at all
    AdaptiveLevel2CenterIndex idx {5, 3};
    auto r = idx.optimistic_bound(10.0, 0, 5.0, false);
    CHECK(r.nterm == 0);
  }

  {
    // Config 6: no discovered, all unknown
    AdaptiveLevel2CenterIndex idx {5, 3};
    auto r = idx.optimistic_bound(10.0, 3, 4.0, false);
    CHECK(r.nterm == 3);
    CHECK(eq(r.density, (10.0 + 3 * 4.0) / 3));
  }
}


TEST_CASE("AdaptiveLevel2ScanCenterIndex matches Fenwick index") {
  using namespace dst;

  AdaptiveLevel2CenterIndex fenwick {5, 6};
  AdaptiveLevel2ScanCenterIndex scan {5, 6};
  std::vector<std::pair<int, double>> discovered {
    {201, 1.0}, {202, 3.0}, {203, 6.0}, {204, 8.0}
  };

  auto compare = [&] (double d_rv, int total, double terminal_lb, bool exhausted) {
    auto f_discovered = fenwick.discovered_bound(d_rv);
    auto s_discovered = scan.discovered_bound(d_rv);
    CHECK(f_discovered.nterm == s_discovered.nterm);
    CHECK(eq(f_discovered.density, s_discovered.density));

    auto f_optimistic = fenwick.optimistic_bound(
        d_rv, total, terminal_lb, exhausted);
    auto s_optimistic = scan.optimistic_bound(
        d_rv, total, terminal_lb, exhausted);
    CHECK(f_optimistic.nterm == s_optimistic.nterm);
    CHECK(eq(f_optimistic.density, s_optimistic.density));
    CHECK(fenwick.discovered_prefix_is_exact(
        d_rv, total, terminal_lb, exhausted) ==
          scan.discovered_prefix_is_exact(
              d_rv, total, terminal_lb, exhausted));

    auto f_tree = fenwick.partree(0, d_rv, f_discovered.nterm);
    auto s_tree = scan.partree(0, d_rv, s_discovered.nterm);
    CHECK(f_tree->terms_before_ready == s_tree->terms_before_ready);
    CHECK(eq(f_tree->cost_sc, s_tree->cost_sc));
  };

  compare(4.0, 5, 2.0, false);
  for (auto [t, d]: discovered) {
    fenwick.insert(t, d);
    scan.insert(t, d);
    compare(4.0, 5, 2.0, false);
  }
  CHECK(scan.best_prefix_length(2.0) == 1); // equal-density extension is omitted
  compare(10.0, 5, 2.0, false); // changing d_rv refreshes both caches

  fenwick.erase({201, 203});
  scan.erase({201, 203});
  compare(4.0, 5, 2.0, false);
  compare(10.0, 2, 0.0, true);

  fenwick.erase({202, 204});
  scan.erase({202, 204});
  compare(10.0, 5, 2.0, false);
}


TEST_CASE("adaptive_fast_level2_no_fenwick matches root-adaptive results") {
  using namespace dst;

  std::vector<std::pair<int,int>> edges {
    {0,1}, {0,2}, {0,3}, {0,4},
    {1,11}, {2,12}, {3,13},
    {4,11}, {4,12}, {4,13}
  };
  std::vector<float> weights {1,1,1,2.5, 2,2,2, 1,1,1};
  DST dt(edges, weights, 0, {11,12,13});

  auto fenwick = dt.adaptive_fast_level2_alg(2);
  auto scan = dt.adaptive_fast_level2_no_fenwick_alg(2);
  auto fenwick_tree = fenwick->to_tree();
  auto scan_tree = scan->to_tree();

  CHECK(eq(fenwick->cost_sc, scan->cost_sc));
  CHECK(eq(fenwick_tree->cost, scan_tree->cost));
  CHECK(fenwick_tree->is_valid(dt.terms_dm));
  CHECK(scan_tree->is_valid(dt.terms_dm));
  CHECK(fenwick->debuginfo.at("center_index") == "fenwick");
  CHECK(scan->debuginfo.at("center_index") == "linear_scan");
  CHECK(scan->debuginfo.at("observed_pool_invariant_violations") == "0");
}


TEST_CASE("adaptive_fast_level2_no_fenwick validates API arguments") {
  using namespace dst;
  std::vector<std::pair<int,int>> edges {{0, 1}, {1, 2}};
  std::vector<float> weights {1, 1};
  DST dt(edges, weights, 0, {2});
  CHECK_THROWS_AS(dt.adaptive_fast_level2_no_fenwick_alg(0), std::invalid_argument);
  CHECK_THROWS_AS(dt.adaptive_fast_level2_no_fenwick_alg(-1), std::invalid_argument);

  DST empty(edges, weights, 0, {});
  auto result = empty.adaptive_fast_level2_no_fenwick_alg();
  CHECK(result->terms.empty());
  CHECK(result->subtrees.empty());
  CHECK(result->debuginfo.at("sssp_nodes_visited") == "0");
  CHECK(result->debuginfo.at("center_index") == "linear_scan");
}


TEST_CASE("full_adaptive_fast_level2 single behavior") {
  using namespace dst;

  /*
              11   12   13
               ^    ^    ^
            (2)|    |    |(2)
               |    |    |
               1    2    3
               ^    ^    ^
            (1)|    |    |(1)
               |    |    |
              <---+0+---+
                    |
                (2.5)
                    v
                    4+--(1)-->11
                    +--(1)-->12
                    +--(1)-->13
  */
  std::vector<std::pair<int,int>> edges {
    {0,1}, {0,2}, {0,3}, {0,4},
    {1,11}, {2,12}, {3,13},
    {4,11}, {4,12}, {4,13}
  };
  std::vector<float> weights {1,1,1,2.5, 2,2,2, 1,1,1};
  std::vector<int> terms {11,12,13};
  DST dt = DST(edges, weights, 0, terms);

  auto faf2 = dt.full_adaptive_fast_level2_alg();

  CHECK(eq(faf2->to_tree()->cost, 5.5));
  CHECK(eq(faf2->cost_sc, 5.5));
}


TEST_CASE("full_adaptive_fast_level2 on adaptive root distance graph") {
  using namespace dst;

  /*
          0
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
  */
  std::vector<std::pair<int,int>> edges {
    {0,1}, {1,11}, {1,2}, {2,21}, {2,22}
  };
  std::vector<float> weights {4,1,5,1,1};
  std::vector<int> terms {11,21,22};
  DST dt = DST(edges, weights, 0, terms);

  auto faf2 = dt.full_adaptive_fast_level2_alg();

  CHECK(eq(faf2->to_tree()->cost, 12));
  CHECK(eq(faf2->cost_sc, 12));
}


TEST_CASE("full_adaptive_fast_level2 on residual-path graph") {
  using namespace dst;

  /*
    0---(1)-->1---(1)-->4---(1)-->6
              +                   ^
              |(2)               (1)
              v                   |
              2---(1)-->5---(1)-->3
              +
              |(3)
              v
              7
  */
  std::vector<std::pair<int,int>> edges {
    {0,1}, {1,4}, {4,6},
    {1,2}, {2,5}, {5,3}, {3,6},
    {2,7}
  };
  std::vector<float> weights {1,1,1, 2,1,1,1, 3};
  std::vector<int> terms {6,7};
  DST dt = DST(edges, weights, 0, terms);

  auto fast2 = dt.fast_level2_alg();
  auto faf2 = dt.full_adaptive_fast_level2_alg();

  // full_adaptive reuses arcs on both sides, so it must do no worse than fast.
  CHECK(leq(faf2->cost_sc, fast2->cost_sc));
  CHECK(eq(faf2->to_tree()->cost, faf2->cost_sc));
}

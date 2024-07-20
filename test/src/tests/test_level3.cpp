#include <fmt/ranges.h>
#include "dst/dst.hpp"
#include "dst/dijkstra.hpp"
#include "dst/utils.hpp"
#include "dst/partree.hpp"

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <tuple>


TEST_CASE("minden_by_thresholds") {
  using namespace dst;

  PartialTreeTable tb1 (1);
  tb1.add_term(11, 1);
  tb1.add_term(12, 10);
  tb1.build();

  PartialTreeTable tb2 (2);
  tb2.add_term(13, 1.5);
  tb2.add_term(14, 1.5);
  tb2.build();

  std::unordered_map<int, PartialTreeTable> tbls;
  tbls[1] = std::move(tb1);
  tbls[2] = std::move(tb2);
  int n_thrs = 3;
  double thr_max = 3;
  auto minden = ThresholdedMinDensity(n_thrs, thr_max, tbls);

  CHECK(minden.thr_mindens->size() == 3);
  CHECK(eq((*minden.thr_mindens)[0].first, 1));
  CHECK(eq((*minden.thr_mindens)[0].second, (1+1)));
  CHECK(eq((*minden.thr_mindens)[1].first, 2));
  CHECK(eq((*minden.thr_mindens)[1].second, (2+1.5*2)/2));
  CHECK(eq((*minden.thr_mindens)[2].first, 3));
  CHECK(eq((*minden.thr_mindens)[2].second, (3+1.5*2)/2));

  CHECK(eq(minden.min_density(0), 0));
  CHECK(eq(minden.min_density(1.5), 2));
  CHECK(eq(minden.min_density(2.5), (2+1.5*2)/2));
  CHECK(eq(minden.min_density(3.5), (3+1.5*2)/2));
}


TEST_CASE("PartialTreeTable") {
  using namespace dst;

  PartialTreeTable tb (10);
  tb.add_term(13, 1.4);
  tb.add_term(11, 1);
  tb.add_term(12, 1);
  tb.add_term(14, 4);
  tb.build();

  CHECK(tb.cells.size() == 4);
  CHECK(eq(tb.cells[0]->d_rv_UB, 0));
  CHECK(tb.find(0.1) == tb.cells[1]);
  CHECK(tb.find(0.79) == tb.cells[1]);
  CHECK(tb.find(0.81) == tb.cells[2]);
  CHECK(tb.find(1) == tb.cells[2]);
  CHECK(tb.find(8.59) == tb.cells[2]);
  CHECK(tb.find(8.61) == tb.cells[3]);
  CHECK(eq(tb.density(0), 1));
  CHECK(eq(tb.density(0.8), 1.4));
  CHECK(eq(tb.density(8.6), 4));
  CHECK(eq(tb.density(10.6), (10.6+1+1+1.4+4)/4));

  PartialTreeTable tb2 (10);
  tb2.add_term(11, 1);
  tb2.add_term(12, 1);
  tb2.add_term(13, 1);
  tb2.add_term(14, 1);
  tb2.build();

  CHECK(tb2.cells.size() == 4);
  CHECK(eq(tb2.cells[0]->d_rv_UB, 0));
  CHECK(tb2.cells[3]->d_rv_UB > 0);
  CHECK(tb2.find(0.1) == tb2.cells[3]);
  CHECK(tb2.find(100) == tb2.cells[3]);

  PartialTreeTable tb3 (10);
  tb3.add_term(11, 1);
  tb3.add_term(12, 2);
  tb3.add_term(13, 3);
  tb3.add_term(14, 4);
  tb3.build();

  CHECK(tb3.cells.size() == 4);
  CHECK(eq(tb3.cells[0]->d_rv_UB, (1.5-1)*2));

  tb3.erase({11,13});
  CHECK(tb3.cells.size() == 2);
  CHECK(eq(tb3.cells[0]->d_rv_UB, (3-2)*2));
}


TEST_CASE("stop_picking_next_2partree") {
  using namespace dst;

  // a full binary tree
  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(0,2), 

                                         std::make_pair(1,11), 
                                         std::make_pair(1,12), 
                                         std::make_pair(2,13), 
                                         std::make_pair(2,14), 

                                         std::make_pair(11,21), 
                                         std::make_pair(11,22), 
                                         std::make_pair(12,23), 
                                         std::make_pair(12,24),
                                         std::make_pair(13,25), 
                                         std::make_pair(13,26), 
                                         std::make_pair(14,27), 
                                         std::make_pair(14,28), 
                                         };
  std::vector<double> weights {1,1, 1,1,1,1, 1,1,1,1,1,1,0.9,100};
  std::vector<int> terms {21,22,23,24,25,26,27,28};
  DST dt = DST(edges, weights, 0, terms);

  auto tree_naive = dt.naive_alg();
  CHECK(tree_naive->cost_sc == 3*6+2.9+2+100);
  CHECK(tree_naive->cost == 2+4+6.9+100);

  auto tree2 = dt.level2_alg(); // pick 11-14 as root 
  CHECK(tree2->to_tree()->cost == 2+4+6.9+100);
  //CHECK(tree2->cost_sc == (1+1+2) * 3 + 2.9 + 2+100);
  CHECK(tree2->cost_sc == (1+1+2) * 3 + 1.9 + 100);

  /*
  next u 1, costsc 7, cov 4
  next u 2, costsc 5.9, cov 3
  next u 14, costsc 102, cov 1
  */
  std::vector<double> costs {7, 5.9, 100};
  std::vector<double> dens {7/4.0, 5.9/3, 100};
  auto tree3 = dt.level3_alg_naive();
  CHECK(tree3->to_tree()->cost == 2+4+6.9+100);
  CHECK(tree3->cost_sc == (1+2+4) + (4+1.9) + 100);
  int i = 0;
  for (auto d: tree3->density1by1) {
    CHECK(eq(d, dens[i])); i++;
  }
  i = 0;
  for (auto c: tree3->costsc1by1) {
    CHECK(eq(c, costs[i])); i++;
  }

  auto fast3 = dt.level3_alg();
  //CHECK(fast3->cost_sc == (1+2+4) + (4+1.9) + 2+100);
  CHECK(fast3->cost_sc == (1+2+4) + (4+1.9) + 100); 
  CHECK(fast3->to_tree()->cost == 2+4+6.9+100);
  i = 0;
  for (auto d: fast3->density1by1) {
    CHECK(eq(d, dens[i])); i++;
  }
  i = 0;
  for (auto c: fast3->costsc1by1) {
    CHECK(eq(c, costs[i])); i++;
  }
}


TEST_CASE("two_3level_partrees") {
  using namespace dst;

  std::vector<std::pair<int,int>> edges {std::make_pair(0,1), 
                                         std::make_pair(0,2), 

                                         std::make_pair(1,11), 
                                         std::make_pair(1,12), 
                                         std::make_pair(2,13), 
                                         std::make_pair(2,14), 

                                         std::make_pair(11,21), 
                                         std::make_pair(11,22), 
                                         std::make_pair(12,23), 
                                         std::make_pair(12,24),
                                         std::make_pair(13,25), 
                                         std::make_pair(13,26), 
                                         std::make_pair(14,27), 
                                         std::make_pair(14,28), 
                                         };
  std::vector<double> weights {1,1, 1,1,1,1, 1,1,1,1,1,1,1,1};
  std::vector<int> terms {21,22,23,24,25,26,27,28};
  DST dt = DST(edges, weights, 0, terms);

  auto tree_naive = dt.naive_alg();
  CHECK(tree_naive->cost_sc == 3*8);
  CHECK(tree_naive->cost == 2+4+8);

  auto tree2 = dt.level2_alg(); // pick 11-14 as root 
  CHECK(tree2->to_tree()->cost == 2+4+8);
  CHECK(tree2->cost_sc == (1+1+2) * 4);

  auto tree3 = dt.level3_alg_naive();
  CHECK(tree3->to_tree()->cost == 2+4+8);
  CHECK(tree3->cost_sc == 2 * (1+2+4));
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

  auto tree_naive = dt.naive_alg();
  CHECK(tree_naive->cost_sc == (2+2) * 4);
  CHECK(tree_naive->cost == (2 + 2*2) * 2);

  auto tree2 = dt.level2_alg();
  CHECK(tree2->to_tree()->cost == (2 + 2*2) * 2);
  CHECK(tree2->cost_sc == (2 + 2*2) * 2);
  CHECK((tree2->terms == std::unordered_set<int> {25,26,27,28}));

  auto tree3 = dt.level3_alg_naive();
  CHECK(tree3->to_tree()->cost == 4 + 2*1.5 + 4*1);
  CHECK(tree3->cost_sc == 4 + 2*1.5 + 4*1);
  CHECK((tree3->terms == std::unordered_set<int> {25,26,27,28}));

  auto fast3 = dt.level3_alg();
  CHECK(fast3->to_tree()->cost == 4 + 2*1.5 + 4*1);
  CHECK(fast3->cost_sc == 4 + 2*1.5 + 4*1);
  CHECK((fast3->terms == std::unordered_set<int> {25,26,27,28}));
}


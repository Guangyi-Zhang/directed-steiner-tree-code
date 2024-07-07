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
}


#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <ctime>

#include "dst/dst.hpp"


int main(int argc, char** argv) {

  using namespace dst;

  /* DATA */

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
  int root = 0;

  /* START RUNNING */

  std::string method = "fast_level2";
  double alpha = 1;
  int nthresholds = 50;

  // one extra run of dijkstra in constructor to remove unreachable vertices
  DST dt = DST(edges, weights, root, terms);

  std::clock_t c_start = std::clock();

  std::shared_ptr<PartialTreeManager> partree = nullptr;
  if (method.compare("level2") == 0) {
    partree = dt.level2_alg();
  }
  else if (method.compare("fast_level2") == 0) {
    partree = dt.fast_level2_alg();
  }
  else if (method.compare("level3") == 0) {
    partree = dt.level3_alg();
  }
  else if (method.compare("fast_level3") == 0) {
    partree = dt.fast_level3_alg(alpha, nthresholds);
  }
  else {
    std::cerr << "unknown method: " << method << std::endl;
    return 1;
  }

  std::clock_t c_end = std::clock();
  double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;

  std::shared_ptr<Tree> tree = partree->to_tree();
  std::cout << "cost " << tree->cost_trimmed() << ", with time " << time_elapsed_ms << " ms" << std::endl;
  tree->print();

  return 0;

}
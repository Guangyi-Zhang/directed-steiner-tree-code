#include <fmt/ranges.h>
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "dst/dst.hpp"

#include <cxxopts.hpp>
#include <cassert>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <tuple>


auto main(int argc, char** argv) -> int {

  using namespace dst;

  // setup logger
  try {
    auto logger = spdlog::basic_logger_mt("basic_logger", "logs/basic-log.txt");
    spdlog::set_default_logger(logger);
  }
  catch (const spdlog::spdlog_ex &ex) {
    std::cout << "Log init failed: " << ex.what() << std::endl;
  }

  // read parameters
  cxxopts::Options options(*argv, "DST");
  std::string buildtype;

  options.add_options()
    ("h,help", "Show help")
    ("b,buildtype", "CMAKE_BUILD_TYPE", cxxopts::value(buildtype)->default_value("unknown"))
  ;
  auto opresult = options.parse(argc, argv);

  // parameters
  std::string version {"v1"}; // @20240616: testing
  int rep {1};

  //std::string method {"naive"};
  std::string method {"level2"};

  double alpha = 0.5;
  std::string dataset {"random_graph_5000.csv"};

  // load data
  std::vector<std::pair<int,int>> edges;
  std::vector<double> weights;

  std::ifstream file("datasets/" + dataset);
  std::string line;
  if (file.is_open()) {
    while (std::getline(file, line)) {
      // split the line into individual values
      std::stringstream ss(line);
      std::string value;

      std::getline(ss, value, ',');
      int v1 = std::stoi(value);
      std::getline(ss, value, ',');
      int v2 = std::stoi(value);
      std::getline(ss, value, ',');
      double w = std::stod(value);

      edges.push_back({v1,v2});
      weights.push_back(w);
    }

    file.close();
  } else {
    std::cerr << "Unable to open file" << std::endl;
    return 1;
  }

  int root {0};
  std::vector<int> terms;
  if (dataset.compare("random_graph_10000.csv") == 0)
    terms = {7270, 860, 5390, 5191, 5734, 6265, 466, 4426, 5578, 8322};
  if (dataset.compare("random_graph_5000.csv") == 0)
    terms = {860, 3772, 3092, 466, 4426, 3444, 3171, 2919, 130, 1685};
  if (dataset.compare("random_graph_1000.csv") == 0)
    terms = {102, 435, 860, 270, 106, 71, 700, 20, 614, 121};
  if (dataset.compare("random_graph_100.csv") == 0)
    terms = {51, 92, 14, 71, 60, 20, 82, 86, 74, 74};
  if (dataset.compare("random_graph_500.csv") == 0)
    terms = {102, 435, 348, 270, 106, 71, 188, 20, 102, 121};
  if (dataset.compare("random_graph_250.csv") == 0)
    terms = {102, 179, 92, 14, 106, 71, 188, 20, 102, 121};
  fmt::println("terms: {}", terms);

  /* START RUNNING */
  std::clock_t c_start = std::clock();

  DST dt = DST(edges, weights, 0, terms);

  double cost {-1}, cost_sc{-1};
  int n_cov {-1};
  if (method.compare("naive") == 0) {
    cost = dt.naive_alg();
  }
  else if (method.compare("level2") == 0) {
    auto tree2 = dt.level2_alg();
    cost = tree2.cost;
    cost_sc = tree2.cost_sc;
    n_cov = tree2.terms_cov.size();
  }
  else if (method.compare("level3") == 0) {
    auto tree3 = dt.level3_alg(alpha);
    cost = tree3.cost;
    cost_sc = tree3.cost_sc;
    n_cov = tree3.terms_cov.size();
  }
  else {
    std::cerr << "unknown method: " << method << std::endl;
    return 1;
  }

  std::clock_t c_end = std::clock();
  double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
  /* FINISH RUNNING */

  // record results as a json string
  rapidjson::Document d;
  rapidjson::Document::AllocatorType& a = d.GetAllocator();

  d.SetObject();
  d.AddMember("version", rapidjson::Value(rapidjson::StringRef(version.c_str())), a); 
  d.AddMember("buildtype", rapidjson::Value(rapidjson::StringRef(opresult["buildtype"].as<std::string>().c_str())), a); 
  d.AddMember("method", rapidjson::Value(rapidjson::StringRef(method.c_str())), a); 
  d.AddMember("rep", rapidjson::Value(rep), a); 
  d.AddMember("dataset", rapidjson::Value(rapidjson::StringRef(dataset.c_str())), a); 
  d.AddMember("alpha", rapidjson::Value(alpha), a); 
  d.AddMember("cost", rapidjson::Value(cost), a); 
  d.AddMember("cost_sc", rapidjson::Value(cost_sc), a); 
  d.AddMember("n_cov", rapidjson::Value(n_cov), a); 
  d.AddMember("runtime", rapidjson::Value(time_elapsed_ms), a); 

  rapidjson::StringBuffer buffer;
  rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
  d.Accept(writer);
  spdlog::info("{}", buffer.GetString());

  return 0;
}

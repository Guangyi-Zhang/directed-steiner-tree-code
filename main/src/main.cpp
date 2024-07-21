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

#include <cxxopts.hpp>
#include <fmt/ranges.h>
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "dst/dst.hpp"
#include "dst/nadeau.hpp"


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

  options.add_options()
    ("h,help", "Show help")
    ("b,buildtype", "CMAKE_BUILD_TYPE", cxxopts::value<std::string>()->default_value("unknown"))
    ("m,method", "method", cxxopts::value<std::string>()->default_value("naive"))
    ("d,dataset", "dataset", cxxopts::value<std::string>()->default_value("random_graph_wFalse_p01_1000.csv"))
    ("v,version", "version", cxxopts::value<int>()->default_value("0"))
    ("r,rep", "rep", cxxopts::value<int>()->default_value("1"))
    ("a,alpha", "alpha", cxxopts::value<double>()->default_value("0.99"))
    ("t,note", "note", cxxopts::value<std::string>()->default_value("none"))
  ;
  auto opresult = options.parse(argc, argv);


  // parameters
  std::string buildtype = opresult["buildtype"].as<std::string>();
  std::string method = opresult["method"].as<std::string>();
  std::string dataset = opresult["dataset"].as<std::string>();
  std::string note = opresult["note"].as<std::string>();
  double alpha = opresult["alpha"].as<double>();
  int rep = opresult["rep"].as<int>();
  int version = opresult["version"].as<int>();

  // load data
  std::vector<std::pair<int,int>> edges;
  std::vector<double> weights;

  std::ifstream file("datasets/" + dataset);
  
  std::string line;
  if (file.is_open()) {
    if (dataset == "soc-Epinions1.txt") {
      while (std::getline(file, line)) {
        if (line[0] == '#')
          continue;

        std::istringstream iss {line};
        std::string value;
        int v1, v2;
        iss >> v1 >> v2;

        edges.push_back({v1,v2});
        weights.push_back(1);
      }
    } else {
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
  if (dataset.find("1000.csv") != std::string::npos)
    terms = {102, 435, 860, 270, 106, 71, 700, 20, 614, 121};
  if (dataset.compare("random_graph_100.csv") == 0)
    terms = {51, 92, 14, 71, 60, 20, 82, 86, 74, 74};
  if (dataset.compare("random_graph_500.csv") == 0)
    terms = {102, 435, 348, 270, 106, 71, 188, 20, 102, 121};
  if (dataset.compare("random_graph_250.csv") == 0)
    terms = {102, 179, 92, 14, 106, 71, 188, 20, 102, 121};
  if (dataset == "soc-Epinions1.txt")
    terms = {7270, 860, 53900, 5191, 25734, 6265, 466, 4426, 65578, 8322};
  fmt::println("terms: {}", terms);

  /* START RUNNING */
  std::clock_t c_start = std::clock();

  DST dt = DST(edges, weights, root, terms);

  std::shared_ptr<Tree> tree = nullptr;
  std::shared_ptr<PartialTreeManager> partree = nullptr;
  std::unordered_map<std::string, std::string> *debuginfo = nullptr;
  if (method.compare("level1") == 0) {
    tree = dt.level1_alg();
    debuginfo = &(tree->debuginfo);
  }
  else if (method.compare("adpnaive_level1") == 0) {
    tree = dt.adaptive_level1_alg();
    debuginfo = &(tree->debuginfo);
  }
  else if (method.compare("level2") == 0) {
    partree = dt.level2_alg();
    tree = partree->to_tree();
    debuginfo = &(partree->debuginfo);
  }
  else if (method.compare("fast_level2") == 0) {
    partree = dt.fast_level2_alg();
    tree = partree->to_tree();
    debuginfo = &(partree->debuginfo);
  }
  else if (method.compare("level3") == 0) {
    partree = dt.level3_alg();
    tree = partree->to_tree();
    debuginfo = &(partree->debuginfo);
  }
  else if (method.compare("fast_level3") == 0) {
    partree = dt.fast_level3_alg(alpha);
    tree = partree->to_tree();
    debuginfo = &(partree->debuginfo);
  }
  else {
    std::cerr << "unknown method: " << method << std::endl;
    return 1;
  }

  std::clock_t c_end = std::clock();
  double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
  /* FINISH RUNNING */

  auto rss = getPeakRSS();

  // record results as a json string
  rapidjson::Document d;
  rapidjson::Document::AllocatorType& a = d.GetAllocator();

  d.SetObject();
  d.AddMember("version", rapidjson::Value(version), a); 
  d.AddMember("note", rapidjson::Value(rapidjson::StringRef(note.c_str())), a); 
  d.AddMember("buildtype", rapidjson::Value(rapidjson::StringRef(opresult["buildtype"].as<std::string>().c_str())), a); 
  d.AddMember("method", rapidjson::Value(rapidjson::StringRef(method.c_str())), a); 
  d.AddMember("rep", rapidjson::Value(rep), a); 
  d.AddMember("dataset", rapidjson::Value(rapidjson::StringRef(dataset.c_str())), a); 
  d.AddMember("alpha", rapidjson::Value(alpha), a); 
  d.AddMember("cost", rapidjson::Value(tree->cost), a); 
  d.AddMember("cost_sc", rapidjson::Value(tree->cost_sc), a); 
  d.AddMember("cost_trimmed", rapidjson::Value(tree->cost_trimmed()), a); 
  d.AddMember("n_cov", rapidjson::Value(tree->terms_cov.size()), a); 
  d.AddMember("runtime", rapidjson::Value(time_elapsed_ms), a); 
  d.AddMember("sssp_nodes_visited", rapidjson::Value(std::stoi(debuginfo->at("sssp_nodes_visited"))), a); 
  d.AddMember("mem", rapidjson::Value(rss), a); 

  rapidjson::StringBuffer buffer;
  rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
  d.Accept(writer);
  spdlog::info("{}", buffer.GetString());

  tree->print();

  return 0;
}

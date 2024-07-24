#include <cassert>
#include <ctime>
#include <cstdlib> // rand()
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

  // read parameters
  cxxopts::Options options(*argv, "DST");

  options.add_options()
    ("h,help", "Show help")
    ("b,buildtype", "CMAKE_BUILD_TYPE", cxxopts::value<std::string>()->default_value("unknown"))
    ("m,method", "method", cxxopts::value<std::string>()->default_value("level1"))
    ("d,dataset", "dataset", cxxopts::value<std::string>()->default_value("none"))
    ("v,version", "version", cxxopts::value<int>()->default_value("0"))
    ("r,rep", "rep", cxxopts::value<int>()->default_value("1"))
    ("s,seed", "seed", cxxopts::value<int>()->default_value("-1"))
    ("k,num_of_terminals", "terminals", cxxopts::value<int>()->default_value("10"))
    ("a,alpha", "alpha", cxxopts::value<double>()->default_value("1"))
    ("t,note", "note", cxxopts::value<std::string>()->default_value("none"))
  ;
  auto opresult = options.parse(argc, argv);


  // parameters
  std::string buildtype = opresult["buildtype"].as<std::string>();
  std::string method = opresult["method"].as<std::string>();
  std::string dataset = opresult["dataset"].as<std::string>();
  std::string note = opresult["note"].as<std::string>();
  double alpha = opresult["alpha"].as<double>();
  int k = opresult["num_of_terminals"].as<int>();
  int rep = opresult["rep"].as<int>();
  int version = opresult["version"].as<int>();
  int seed = opresult["seed"].as<int>();
  if (seed == -1)
    seed = time(nullptr);
  srand(seed);
  spdlog::info("running method={}, dataset={}, ver={}, rep={}, seed={}", method, dataset, version, rep, seed);

  // load data
  std::vector<std::pair<int,int>> edges;
  std::vector<double> weights;

  std::ifstream file("datasets/" + dataset);
  
  std::string line;
  int n = 0, m = 0;
  std::unordered_set<int> V; // vertex id in [0, V.size())
  if (file.is_open()) {
    if (dataset == "soc-Epinions1.txt" or 
        dataset == "web-Google.txt" or
        dataset == "soc-pokec-relationships.txt"
       ) 
    {
      std::unordered_map<int,int> id2v;
      while (std::getline(file, line)) {
        if (line[0] == '#')
          continue;

        std::istringstream iss {line};
        int v1, v2;
        iss >> v1 >> v2;

        if (not has_key(id2v, v1))
          id2v[v1] = id2v.size();
        v1 = id2v.at(v1);
        if (not has_key(id2v, v2))
          id2v[v2] = id2v.size();
        v2 = id2v.at(v2);

        edges.push_back({v1,v2});
        weights.push_back(1);
        V.insert(v1); V.insert(v2);
        m++;
      }
    } 
    else if (dataset == "SFRoad") {
      std::unordered_map<int,int> id2v;
      while (std::getline(file, line)) {
        std::istringstream iss {line};
        int edgeid, v1, v2;
        double w;
        iss >> edgeid >> v1 >> v2 >> w;

        if (not has_key(id2v, v1))
          id2v[v1] = id2v.size();
        v1 = id2v.at(v1);
        if (not has_key(id2v, v2))
          id2v[v2] = id2v.size();
        v2 = id2v.at(v2);

        edges.push_back({v1,v2});
        weights.push_back(w);
        V.insert(v1); V.insert(v2);
        m++;
      }
    }
    else if (dataset == "token_transfers.csv") {
      std::unordered_map<std::string,int> addr2v;
      std::getline(file, line); // skip 1st line
      while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        std::getline(ss, value, ','); // block_number
        std::getline(ss, value, ','); // transaction_index
        std::getline(ss, value, ','); // from_address
        if (not has_key(addr2v, value))
          addr2v[value] = addr2v.size();
        int v1 = addr2v.at(value);
        std::getline(ss, value, ','); // to_address
        if (not has_key(addr2v, value))
          addr2v[value] = addr2v.size();
        int v2 = addr2v.at(value);
        std::getline(ss, value, ','); // time_stamp
        std::getline(ss, value, ','); // contract_address
        std::getline(ss, value, ','); // value
        double w = std::stod(value);

        edges.push_back({v1,v2});
        weights.push_back(w);
        V.insert(v1); V.insert(v2);
        m++;
      }
    }
    else { // synthetic random data
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
        V.insert(v1); V.insert(v2);
        m++;
      }
    }
    n = V.size();

    file.close();
  } else {
    spdlog::error("Error: unable to open file");
    return 1;
  }

  int root = rand() % n;
  std::vector<int> terms;
  std::unordered_set<int> terms_set;
  spdlog::info("RAND_MAX: {}", RAND_MAX); // same as INT_MAX 
  for (int i=0; i<k; i++) {
    // generate random terminals, [0, n)
    int t = rand() % n;
    while(has_key(terms_set, t))
      t = rand() % n;
    terms.push_back(t);
    terms_set.insert(t);
  }
  spdlog::info("root={}, terms= {}", root, terms);


  /* START RUNNING */

  // one extra run of dijkstra in constructor to remove unreachable vertices
  DST dt = DST(edges, weights, root, terms);

  std::clock_t c_start = std::clock();

  std::shared_ptr<Tree> tree = nullptr;
  std::shared_ptr<PartialTreeManager> partree = nullptr;
  std::unordered_map<std::string, std::string> *debuginfo = nullptr;
  if (method.compare("level1") == 0) {
    tree = dt.level1_alg();
    debuginfo = &(tree->debuginfo);
  }
  else if (method.compare("adaptive_level1") == 0) {
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
  d.AddMember("buildtype", rapidjson::Value(rapidjson::StringRef(buildtype.c_str())), a); 
  d.AddMember("method", rapidjson::Value(rapidjson::StringRef(method.c_str())), a); 
  d.AddMember("rep", rapidjson::Value(rep), a); 
  d.AddMember("k", rapidjson::Value(k), a); 
  d.AddMember("n", rapidjson::Value(n), a); 
  d.AddMember("m", rapidjson::Value(m), a); 
  d.AddMember("k_", rapidjson::Value(dt.terms.size()), a); 
  d.AddMember("n_", rapidjson::Value(dt.V.size()), a); 
  d.AddMember("m_", rapidjson::Value(dt.narcs), a); 
  d.AddMember("seed", rapidjson::Value(seed), a); 
  d.AddMember("root", rapidjson::Value(root), a); 
  d.AddMember("dataset", rapidjson::Value(rapidjson::StringRef(dataset.c_str())), a); 
  d.AddMember("alpha", rapidjson::Value(alpha), a); 
  d.AddMember("cost", rapidjson::Value(tree->cost), a); 
  d.AddMember("cost_sc", rapidjson::Value(tree->cost_sc), a); 
  d.AddMember("cost_trimmed", rapidjson::Value(tree->cost_trimmed()), a); 
  d.AddMember("n_cov", rapidjson::Value(tree->terms_cov.size()), a); 
  d.AddMember("runtime", rapidjson::Value(time_elapsed_ms), a); 
  d.AddMember("sssp_nodes_visited", rapidjson::Value(std::stoi(debuginfo->at("sssp_nodes_visited"))), a); 
  d.AddMember("mem", rapidjson::Value(rss), a); 

  std::stringstream ss_terms;
  std::copy(terms.begin(), terms.end(), std::ostream_iterator<int>(ss_terms, ","));
  std::string sterms = ss_terms.str();
  sterms.erase(sterms.length()-1);
  rapidjson::Value vterms;
  vterms.SetString(sterms.c_str(), sterms.length(), a); // rapidjson::UTF8 not compatible with string
  d.AddMember("terms", vterms, a); 

  rapidjson::StringBuffer buffer;
  rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
  d.Accept(writer);

  // setup logger
  try {
    auto logger = spdlog::basic_logger_mt("basic_logger", "logs/basic-log.txt");
    spdlog::set_default_logger(logger);
  }
  catch (const spdlog::spdlog_ex &ex) {
    spdlog::error("Log init failed: {}", ex.what());
  }
  spdlog::info("{}", buffer.GetString()); // re-directed to file

  tree->print();

  return 0;
}

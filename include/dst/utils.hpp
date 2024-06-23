#pragma once

#include <numeric> // std::iota
#include <algorithm>
#include <iterator> // inserter
#include <vector>
#include <functional> // std::hash


namespace dst {
  struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
      auto h1 = std::hash<T1>{}(p.first);
      auto h2 = std::hash<T2>{}(p.second);

      // Mainly for demonstration purposes, i.e. works but is overly simple
      // In the real world, use sth. like boost.hash_combine
      return h1 ^ h2;  
    }
  };
  

  template <typename T, typename U>
  bool has_key(const T& container, const U& key) {
    // only works for set and map!
    return (container.find(key) != container.end());
  }


  template <typename T>
  std::vector<size_t> argsort(const std::vector<T> &v) {
    // https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0); // fill idx by increasing values

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    stable_sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
  }

  
}
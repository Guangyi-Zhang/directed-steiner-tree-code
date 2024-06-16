#pragma once

#include <numeric> // std::iota
#include <algorithm>
#include <iterator> // inserter
#include <vector>


namespace dst {

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
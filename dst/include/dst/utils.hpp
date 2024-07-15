#pragma once

#include <cmath>
#include <numeric> // std::iota
#include <algorithm>
#include <iterator> // inserter
#include <vector>
#include <functional> // std::hash

#include <dst/consts.hpp>


namespace dst {

  bool eq(double a, double b) {
    return std::abs(a-b) < EPSILON;
  }

  bool lq(double a, double b) {
    if (std::abs(a-b) < EPSILON)
      return false;
    return a < b;
  }

  bool leq(double a, double b) {
    return a<b or std::abs(a-b) < EPSILON;
  }

  template <typename T, typename U>
  bool has_key(const T& container, const U& key) {
    // only works for set and map!
    return (container.find(key) != container.end());
  }


  // a hack to iterate priority_queue: 
  // https://stackoverflow.com/questions/4484767/how-to-iterate-over-a-priority-queue
  template <class T, class S, class C>
  S& Container(std::priority_queue<T, S, C>& q) {
    struct HackedQueue : private std::priority_queue<T, S, C> {
      static S& Container(std::priority_queue<T, S, C>& q) {
        return q.*&HackedQueue::c;
      }
    };

    return HackedQueue::Container(q);
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
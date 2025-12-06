#pragma once
#include <chrono>

template <typename T>
long benchmark_time(T func, int iterations=25) {
  std::vector<double> vec_diff;
  vec_diff.reserve(iterations);
  for (int i = 0; i < iterations; ++i) {
    auto start = std::chrono::steady_clock::now();
    func();
    auto end = std::chrono::steady_clock::now();

    // Convert to nanoseconds
    double diff = std::chrono::duration<double, std::nano>(end - start).count();
    vec_diff.push_back(diff);
    }

  // get average
  double avg = 0.0;
  for (double val : vec_diff) {
      avg += val;
  }
  // return time diff in microseconds
  return avg / vec_diff.size();
}

// need to implement L2 norm for this benchmarking
template <typename T>
long benchmark_acc(T func1, T func2) {
  auto diff = func1() - func2();
  return diff;
}
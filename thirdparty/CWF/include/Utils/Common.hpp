#pragma once
#include <chrono>
#include <spdlog/spdlog.h>
#include <string>

namespace Utils {
class ScopeTimer {
public:
  std::string _name;
  std::chrono::high_resolution_clock::time_point _start;
  bool stopped = false;

  ScopeTimer(const std::string &name) : _name(name) { _start = std::chrono::high_resolution_clock::now(); }

  ~ScopeTimer() { Print(); }

  void Print() {
    if (stopped)
      return;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - _start);
    spdlog::info("耗时: {:2.4f}s [{}]", duration.count() * 0.001, _name);
  }

  void Stop() {
    Print();
    stopped = true;
  }
};
} // namespace Utils
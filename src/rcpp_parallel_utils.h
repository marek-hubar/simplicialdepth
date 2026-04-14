#pragma once

#include <RcppParallel.h>

#include <cstddef>
#include <utility>

// Optional (per-call) thread control.
// threads <= 0 means: use RcppParallel/TBB default.
#if __has_include(<tbb/task_arena.h>)
  #include <tbb/task_arena.h>
  #define SIMPLICIALDEPTH_HAS_TBB_TASK_ARENA 1
#else
  #define SIMPLICIALDEPTH_HAS_TBB_TASK_ARENA 0
#endif

namespace simplicialdepth {

inline int normalize_threads(int threads) {
    return threads <= 0 ? 0 : threads;
}

template <typename Func>
inline void with_threads(int threads, Func&& f) {
#if SIMPLICIALDEPTH_HAS_TBB_TASK_ARENA
    const int normalized = normalize_threads(threads);
    if (normalized > 0) {
        tbb::task_arena arena(normalized);
        arena.execute(std::forward<Func>(f));
        return;
    }
#endif
    std::forward<Func>(f)();
}

// Small wrapper to keep call sites clean.
template <typename Body>
inline void parallel_for(std::size_t begin, std::size_t end, Body& body, int threads = 0) {
    with_threads(threads, [&]() {
        RcppParallel::parallelFor(static_cast<std::size_t>(begin), static_cast<std::size_t>(end), body);
    });
}

template <typename Reducer>
inline void parallel_reduce(std::size_t begin,
                            std::size_t end,
                            Reducer& reducer,
                            std::size_t grainSize = 1,
                            int threads = 0) {
    with_threads(threads, [&]() {
        RcppParallel::parallelReduce(begin, end, reducer, grainSize);
    });
}

} // namespace simplicialdepth

#pragma once

#include <vector>

namespace simplicialdepth {

struct SweepLine {
    double m;
    double b;
};

struct SweepIntersectionPair {
    int first;
    int second;

    SweepIntersectionPair(int a, int b);
};

class TopologicalSweep {
public:
    std::vector<SweepIntersectionPair> sweep(const std::vector<SweepLine>& lines) const;
};

std::vector<int> build_order_matrix(const std::vector<SweepIntersectionPair>& intersections,
                                    int line_count);

} // namespace simplicialdepth


#include "topological_sweep.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace simplicialdepth {

SweepIntersectionPair::SweepIntersectionPair(int a, int b)
    : first(std::min(a, b)), second(std::max(a, b)) {}

namespace {

struct HTLine {
    int line_index;
    int upper_parent_index;
    int lower_parent_index;
};

double intersection_x(const SweepLine& line1, const SweepLine& line2) {
    return (line1.b - line2.b) / (line2.m - line1.m);
}

class Order {
public:
    void init(int n) {
        order.resize(static_cast<std::size_t>(n));
        ranks_of_ids.resize(static_cast<std::size_t>(n));
        std::iota(order.begin(), order.end(), 0);
        std::iota(ranks_of_ids.begin(), ranks_of_ids.end(), 0);
    }

    int at(int position) const {
        return order[static_cast<std::size_t>(position)];
    }

    int rank(int id) const {
        return ranks_of_ids[static_cast<std::size_t>(id)];
    }

    void swap_ids(int id1, int id2) {
        const int rank1 = ranks_of_ids[static_cast<std::size_t>(id1)];
        const int rank2 = ranks_of_ids[static_cast<std::size_t>(id2)];
        std::swap(order[static_cast<std::size_t>(rank1)], order[static_cast<std::size_t>(rank2)]);
        std::swap(ranks_of_ids[static_cast<std::size_t>(id1)], ranks_of_ids[static_cast<std::size_t>(id2)]);
    }

private:
    std::vector<int> order;
    std::vector<int> ranks_of_ids;
};

void init_horizon_trees(const std::vector<SweepLine>& lines,
                        std::vector<HTLine>& htlines,
                        Order& order) {
    const int n = static_cast<int>(lines.size());
    htlines.resize(static_cast<std::size_t>(n));
    order.init(n);

    for (int i = 0; i < n; i++) {
        htlines[static_cast<std::size_t>(i)] = {i, -1, -1};
    }

    std::sort(htlines.begin(), htlines.end(), [&lines](const HTLine& a, const HTLine& b) {
        const SweepLine& la = lines[static_cast<std::size_t>(a.line_index)];
        const SweepLine& lb = lines[static_cast<std::size_t>(b.line_index)];
        if (la.m != lb.m) return la.m >= lb.m;
        return la.b >= lb.b;
    });

    for (int i = 1; i < n; i++) {
        int j = i - 1;
        double prev_x_intercept = intersection_x(
            lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(j)].line_index)],
            lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(i)].line_index)]);
        int prev_j = j;
        j = htlines[static_cast<std::size_t>(j)].upper_parent_index;
        while (j >= 0) {
            const double curr_x_intercept = intersection_x(
                lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(j)].line_index)],
                lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(i)].line_index)]);
            if (curr_x_intercept > prev_x_intercept) break;
            prev_x_intercept = curr_x_intercept;
            prev_j = j;
            j = htlines[static_cast<std::size_t>(j)].upper_parent_index;
        }

        htlines[static_cast<std::size_t>(i)].upper_parent_index = prev_j;
    }

    for (int i = n - 2; i >= 0; i--) {
        int j = i + 1;
        double prev_x_intercept = intersection_x(
            lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(j)].line_index)],
            lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(i)].line_index)]);
        int prev_j = j;
        j = htlines[static_cast<std::size_t>(j)].lower_parent_index;
        while (j >= 0 && j < n) {
            const double curr_x_intercept = intersection_x(
                lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(j)].line_index)],
                lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(i)].line_index)]);
            if (curr_x_intercept > prev_x_intercept) break;
            prev_x_intercept = curr_x_intercept;
            prev_j = j;
            j = htlines[static_cast<std::size_t>(j)].lower_parent_index;
        }

        htlines[static_cast<std::size_t>(i)].lower_parent_index = prev_j;
    }
}

} // namespace

std::vector<SweepIntersectionPair> TopologicalSweep::sweep(const std::vector<SweepLine>& lines) const {
    const int n = static_cast<int>(lines.size());
    if (n < 2) return {};

    std::vector<SweepIntersectionPair> intersection_order;
    intersection_order.reserve(static_cast<std::size_t>(n * (n - 1) / 2));

    std::vector<HTLine> htlines;
    Order order;
    init_horizon_trees(lines, htlines, order);

    std::vector<int> stack;
    stack.reserve(static_cast<std::size_t>(n));
    for (int i = 1; i < n; i++) {
        if (htlines[static_cast<std::size_t>(i)].upper_parent_index == i - 1 &&
            htlines[static_cast<std::size_t>(i - 1)].lower_parent_index == i) {
            stack.push_back(i);
        }
    }

    while (!stack.empty()) {
        const int popped_idx = stack.back();
        stack.pop_back();

        const HTLine popped = htlines[static_cast<std::size_t>(popped_idx)];
        const int higher_slope_index = popped.upper_parent_index;
        const int lower_slope_index = popped_idx;
        if (higher_slope_index < 0 || higher_slope_index >= n) continue;

        intersection_order.emplace_back(
            htlines[static_cast<std::size_t>(lower_slope_index)].line_index,
            htlines[static_cast<std::size_t>(higher_slope_index)].line_index);

        order.swap_ids(lower_slope_index, higher_slope_index);

        const double x_coordinate_of_intersection = intersection_x(
            lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(higher_slope_index)].line_index)],
            lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(lower_slope_index)].line_index)]);

        int rank_lower = order.rank(lower_slope_index);
        int rank_higher = order.rank(higher_slope_index);

        if (rank_lower > 0) {
            int i = order.at(rank_lower - 1);
            double curr_x_intercept = intersection_x(
                lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(i)].line_index)],
                lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(lower_slope_index)].line_index)]);
            double prev_x_intercept = -std::numeric_limits<double>::infinity();
            int prev_i = i;

            while (!(prev_x_intercept < curr_x_intercept &&
                     x_coordinate_of_intersection < prev_x_intercept)) {
                prev_i = i;
                i = htlines[static_cast<std::size_t>(i)].upper_parent_index;
                if (i == -1) break;
                prev_x_intercept = curr_x_intercept;
                curr_x_intercept = intersection_x(
                    lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(i)].line_index)],
                    lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(lower_slope_index)].line_index)]);
            }

            if (i == -1 &&
                intersection_x(
                    lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(prev_i)].line_index)],
                    lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(lower_slope_index)].line_index)]) <
                    x_coordinate_of_intersection) {
                htlines[static_cast<std::size_t>(lower_slope_index)].upper_parent_index = -1;
            } else {
                htlines[static_cast<std::size_t>(lower_slope_index)].upper_parent_index = prev_i;
            }
        } else {
            htlines[static_cast<std::size_t>(lower_slope_index)].upper_parent_index = -1;
        }

        if (rank_higher < n - 1) {
            int i = order.at(rank_higher + 1);
            double curr_x_intercept = intersection_x(
                lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(i)].line_index)],
                lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(higher_slope_index)].line_index)]);
            double prev_x_intercept = -std::numeric_limits<double>::infinity();
            int prev_i = i;

            while (!(prev_x_intercept < curr_x_intercept &&
                     x_coordinate_of_intersection < prev_x_intercept)) {
                prev_i = i;
                i = htlines[static_cast<std::size_t>(i)].lower_parent_index;
                if (i == -1) break;
                prev_x_intercept = curr_x_intercept;
                curr_x_intercept = intersection_x(
                    lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(i)].line_index)],
                    lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(higher_slope_index)].line_index)]);
            }

            if (i == -1 &&
                intersection_x(
                    lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(prev_i)].line_index)],
                    lines[static_cast<std::size_t>(htlines[static_cast<std::size_t>(higher_slope_index)].line_index)]) <
                    x_coordinate_of_intersection) {
                htlines[static_cast<std::size_t>(higher_slope_index)].lower_parent_index = -1;
            } else {
                htlines[static_cast<std::size_t>(higher_slope_index)].lower_parent_index = prev_i;
            }
        } else {
            htlines[static_cast<std::size_t>(higher_slope_index)].lower_parent_index = -1;
        }

        if (rank_higher + 1 < n) {
            const int index_above_higher_slope = order.at(rank_higher + 1);
            if (htlines[static_cast<std::size_t>(index_above_higher_slope)].upper_parent_index ==
                    higher_slope_index &&
                htlines[static_cast<std::size_t>(higher_slope_index)].lower_parent_index ==
                    index_above_higher_slope) {
                stack.push_back(index_above_higher_slope);
            }
        }

        if (rank_lower - 1 >= 0) {
            const int index_below_lower_slope = order.at(rank_lower - 1);
            if (htlines[static_cast<std::size_t>(lower_slope_index)].upper_parent_index ==
                    index_below_lower_slope &&
                htlines[static_cast<std::size_t>(index_below_lower_slope)].lower_parent_index ==
                    lower_slope_index) {
                stack.push_back(lower_slope_index);
            }
        }
    }

    return intersection_order;
}

std::vector<int> build_order_matrix(const std::vector<SweepIntersectionPair>& intersections,
                                    int line_count) {
    if (line_count < 0) {
        throw std::invalid_argument("line_count must be non-negative");
    }
    if (line_count < 2) return {};

    const std::size_t expected_intersections =
        static_cast<std::size_t>(line_count) * static_cast<std::size_t>(line_count - 1) / 2;
    if (intersections.size() != expected_intersections) {
        throw std::invalid_argument("intersection count does not match line_count");
    }

    const int row_width = line_count - 1;
    std::vector<int> order_matrix(static_cast<std::size_t>(line_count * row_width), -1);
    std::vector<int> indices(static_cast<std::size_t>(line_count), 0);

    for (const auto& pair : intersections) {
        if (pair.first < 0 || pair.first >= line_count || pair.second < 0 || pair.second >= line_count ||
            pair.first == pair.second) {
            throw std::invalid_argument("intersection pair contains invalid line index");
        }

        const int first_col = indices[static_cast<std::size_t>(pair.first)]++;
        const int second_col = indices[static_cast<std::size_t>(pair.second)]++;
        if (first_col >= row_width || second_col >= row_width) {
            throw std::runtime_error("intersection pair overflow while building order matrix");
        }

        order_matrix[static_cast<std::size_t>(pair.first * row_width + first_col)] = pair.second;
        order_matrix[static_cast<std::size_t>(pair.second * row_width + second_col)] = pair.first;
    }

    for (int i = 0; i < line_count; i++) {
        if (indices[static_cast<std::size_t>(i)] != row_width) {
            throw std::runtime_error("incomplete order matrix construction");
        }
    }

    return order_matrix;
}

} // namespace simplicialdepth


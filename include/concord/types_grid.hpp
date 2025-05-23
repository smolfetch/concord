
#pragma once

#include "types_basic.hpp"
#include "types_polygon.hpp"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace concord {

    template <typename T> class Grid {
        using size_type = std::size_t;
        using value_type = std::pair<Point, T>;
        using const_reference = const value_type &;
        using reference = value_type &;

      private:
        size_type rows_{0};
        size_type cols_{0};
        double inradius_{0.0};
        std::vector<value_type> data_;

      public:
        Grid() = default;

        Grid(size_type rows, size_type cols, double diameter, concord::Datum datum, bool centered, concord::Pose shift)
            : rows_{rows}, cols_{cols}, inradius_{diameter}, data_(rows * cols) {
            // Total span of the grid (only non-zero if centered)
            double height = centered ? cols_ * diameter : 0.0;
            double width = centered ? rows_ * diameter : 0.0;
            // Precompute rotation sin/cos
            double yaw = shift.angle.yaw;
            double cosYaw = std::cos(yaw);
            double sinYaw = std::sin(yaw);
            for (size_type r = 0; r < rows_; ++r) {
                for (size_type c = 0; c < cols_; ++c) {
                    // Compute each cell center in the unrotated, unshifted local frame
                    double local_x = ((static_cast<double>(r) + 0.5) * diameter) - (width / 2.0);
                    double local_y = ((static_cast<double>(c) + 0.5) * diameter) - (height / 2.0);
                    // Rotate around origin by yaw
                    double rot_x = cosYaw * local_x - sinYaw * local_y;
                    double rot_y = sinYaw * local_x + cosYaw * local_y;
                    // Then translate by the shift’s ENU offset
                    double world_x = rot_x + shift.point.enu.x;
                    double world_y = rot_y + shift.point.enu.y;
                    // Fill in the point and convert to WGS
                    size_type idx = index(r, c);
                    Point p;
                    p.enu.x = world_x;
                    p.enu.y = world_y;
                    p.enu.toWGS(datum);
                    data_[idx] = value_type{p, T{}};
                }
            }
        }

        reference operator()(size_type r, size_type c) noexcept { return data_[index(r, c)]; }
        const_reference operator()(size_type r, size_type c) const noexcept { return data_[index(r, c)]; }
        reference at(size_type r, size_type c) { return data_.at(index_checked(r, c)); }
        const_reference at(size_type r, size_type c) const { return data_.at(index_checked(r, c)); }
        std::span<value_type> row(size_type r) noexcept { return {&data_[index(r, 0)], cols_}; }
        std::span<const value_type> row(size_type r) const noexcept { return {&data_[index(r, 0)], cols_}; }
        constexpr size_type index(size_type r, size_type c) const noexcept { return r * cols_ + c; }

        size_type index_checked(size_type r, size_type c) const {
            if (r >= rows_ || c >= cols_) {
                throw std::out_of_range("Grid indices out of bounds");
            }
            return r * cols_ + c;
        }

        constexpr size_type rows() const noexcept { return rows_; }
        constexpr size_type cols() const noexcept { return cols_; }
        constexpr double inradius() const noexcept { return inradius_; }

        auto begin() noexcept { return data_.begin(); }
        auto end() noexcept { return data_.end(); }
        auto begin() const noexcept { return data_.begin(); }
        auto end() const noexcept { return data_.end(); }

        std::vector<std::array<float, 3>> flatten_points() {
            std::vector<std::array<float, 3>> points;
            for (auto &[p, c] : data_) {
                points.push_back({float(p.enu.x), float(p.enu.y), 0.0f});
            }
            return points;
        }

        std::vector<std::array<float, 3>> flatten_points(bool just_bool) {
            std::vector<std::array<float, 3>> points;
            for (std::size_t r = 0; r < rows_; ++r) {
                for (std::size_t c = 0; c < cols_; ++c) {
                    auto &[pt, color] = data_[index(r, c)]; // now color is an RGB& directly
                    points.push_back({float(pt.enu.x), float(pt.enu.y), 0.0f});
                }
            }
            return points;
        }

        std::array<concord::Point, 4> corners(concord::Datum datum = Datum()) const {
            if (rows_ == 0 || cols_ == 0) {
                throw std::runtime_error("Grid is empty; cannot get corners");
            }
            const std::size_t r0 = 0, r1 = rows_ - 1;
            const std::size_t c0 = 0, c1 = cols_ - 1;

            auto getP = [&](std::size_t r, std::size_t c) {
                Point p = (*this)(r, c).first;
                p.wgs = p.enu.toWGS(datum);
                return p;
            };

            return {
                getP(r0, c0), // top-left
                getP(r0, c1), // top-right
                getP(r1, c1), // bottom-right
                getP(r1, c0)  // bottom-left
            };
        }

        std::vector<size_type> indices_within(const Polygon &poly) const {
            if (!poly.isConnected()) {
                throw std::invalid_argument("Polygon must have at least 3 vertices");
            }
            std::vector<size_type> out;
            out.reserve(rows_ * cols_);
            for (size_type r = 0; r < rows_; ++r) {
                for (size_type c = 0; c < cols_; ++c) {
                    size_type idx = index(r, c);
                    if (poly.contains(data_[idx].first)) {
                        out.push_back(idx);
                    }
                }
            }
            return out;
        }
    };

} // namespace concord

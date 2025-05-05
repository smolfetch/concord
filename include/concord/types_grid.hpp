#pragma once

#include "types_basic.hpp"
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

        Grid(size_type rows, size_type cols, double inradius)
            : rows_{rows}, cols_{cols}, inradius_{inradius}, data_(rows * cols) {
            double diameter = 2.0 * inradius_;
            for (size_type r = 0; r < rows_; ++r) {
                for (size_type c = 0; c < cols_; ++c) {
                    double center_x = (static_cast<double>(c) + 0.5) * diameter;
                    double center_y = (static_cast<double>(r) + 0.5) * diameter;
                    size_type idx = index(r, c);
                    Point p;
                    p.enu.x = center_x;
                    p.enu.y = center_y;
                    data_[idx] = value_type{p, T{}};
                }
            }
        }

        Grid(size_type rows, size_type cols, double inradius, concord::Datum datum)
            : rows_{rows}, cols_{cols}, inradius_{inradius}, data_(rows * cols) {
            double diameter = 2.0 * inradius_;
            for (size_type r = 0; r < rows_; ++r) {
                for (size_type c = 0; c < cols_; ++c) {
                    double center_x = (static_cast<double>(c) + 0.5) * diameter;
                    double center_y = (static_cast<double>(r) + 0.5) * diameter;
                    size_type idx = index(r, c);
                    Point p;
                    p.enu.x = center_x;
                    p.enu.y = center_y;
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

        constexpr size_type rows() const noexcept { return rows_; }
        constexpr size_type cols() const noexcept { return cols_; }
        constexpr double inradius() const noexcept { return inradius_; }

        auto begin() noexcept { return data_.begin(); }
        auto end() noexcept { return data_.end(); }
        auto begin() const noexcept { return data_.begin(); }
        auto end() const noexcept { return data_.end(); }

      private:
        constexpr size_type index(size_type r, size_type c) const noexcept { return r * cols_ + c; }

        size_type index_checked(size_type r, size_type c) const {
            if (r >= rows_ || c >= cols_) {
                throw std::out_of_range("Grid indices out of bounds");
            }
            return r * cols_ + c;
        }
    };

} // namespace concord

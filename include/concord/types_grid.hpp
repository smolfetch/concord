#pragma once

#include "types_basic.hpp"
#include <cassert>
#include <cstddef>
#include <span> // C++20
#include <vector>

namespace concord {

    class Grid {
      public:
        using size_type = std::size_t;

        // constructor: allocate and fill in one shot
        Grid(size_type rows, size_type cols, double inradius)
            : rows_{rows}, cols_{cols}, inradius_{inradius},
              data_{rows * cols, Point{ENU(inradius, inradius, 0.0), WGS{0.0, 0.0, 0.0}}} {}

        // unchecked access
        Point &operator()(size_type r, size_type c) noexcept { return data_[index(r, c)]; }
        Point const &operator()(size_type r, size_type c) const noexcept { return data_[index(r, c)]; }

        // bounds‐checking version
        Point &at(size_type r, size_type c) { return data_.at(index_checked(r, c)); }
        Point const &at(size_type r, size_type c) const { return data_.at(index_checked(r, c)); }

        // give me a full row as a span
        std::span<Point> row(size_type r) noexcept { return {&data_[index(r, 0)], cols_}; }
        std::span<Point const> row(size_type r) const noexcept { return {&data_[index(r, 0)], cols_}; }

        // simple getters
        constexpr size_type rows() const noexcept { return rows_; }
        constexpr size_type cols() const noexcept { return cols_; }
        constexpr double inradius() const noexcept { return inradius_; }

        // iteration
        auto begin() noexcept { return data_.begin(); }
        auto end() noexcept { return data_.end(); }
        auto begin() const noexcept { return data_.begin(); }
        auto end() const noexcept { return data_.end(); }

      private:
        // plain index into the flat array
        constexpr size_type index(size_type r, size_type c) const noexcept { return r * cols_ + c; }

        // same, but assert in debug builds
        size_type index_checked(size_type r, size_type c) const {
            assert(r < rows_ && c < cols_);
            return r * cols_ + c;
        }

        size_type rows_, cols_;
        double inradius_;
        std::vector<Point> data_;
    };

} // namespace concord

#pragma once
#include "types_basic.hpp"
#include <cassert> // for assert()
#include <cstddef> // for std::size_t
#include <vector>

namespace concord {

    class Grid {
      public:
        // constructor: set size and inradius, and allocate storage
        Grid(std::size_t rows, std::size_t cols, double inradius)
            : rows_{rows}, cols_{cols}, inradius_{inradius}, data_(rows * cols) {}

        // element access, no bounds-checking
        Point &at(std::size_t r, std::size_t c) { return data_[index(r, c)]; }
        const Point &at(std::size_t r, std::size_t c) const { return data_[index(r, c)]; }

        // operator() for nicer syntax
        Point &operator()(std::size_t r, std::size_t c) { return at(r, c); }
        const Point &operator()(std::size_t r, std::size_t c) const { return at(r, c); }

        // optionally, a checked access
        Point &safe_at(std::size_t r, std::size_t c) {
            assert(r < rows_ && c < cols_);
            return data_[r * cols_ + c];
        }

        // getters
        std::size_t rows() const { return rows_; }
        std::size_t cols() const { return cols_; }
        double inradius() const { return inradius_; }

        // iteration support
        auto begin() { return data_.begin(); }
        auto end() { return data_.end(); }
        auto begin() const { return data_.begin(); }
        auto end() const { return data_.end(); }

      private:
        std::size_t index(std::size_t r, std::size_t c) const { return r * cols_ + c; }

        std::size_t rows_, cols_;
        double inradius_;
        std::vector<Point> data_;
    };
} // namespace concord

#pragma once

#include "types_basic.hpp"
#include <cmath>

namespace concord {

    class Square {
      public:
        Square() = default;
        Square(const Point &c, double s) : center(c), side(s) {}

        double area() const noexcept { return side * side; }
        double perimeter() const noexcept { return 4 * side; }
        double diagonal() const noexcept { return side * std::sqrt(2.0); }

        bool contains(const Point &p) const noexcept {
            return std::abs(p.enu.x - center.enu.x) <= side / 2.0 && std::abs(p.enu.y - center.enu.y) <= side / 2.0;
        }

        const Point &getCenter() const noexcept { return center; }
        double getSide() const noexcept { return side; }

      private:
        Point center;
        double side;
    };

} // namespace concord

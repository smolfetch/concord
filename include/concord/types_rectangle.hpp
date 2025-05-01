#pragma once

#include "types_basic.hpp"
#include "types_line.hpp"
#include <cstddef>

namespace concord {

    class Rectangle {
      public:
        Rectangle(const Point &tl, const Point &tr, const Point &bl, const Point &br)
            : top_left(tl), top_right(tr), bottom_left(bl), bottom_right(br) {}

        double area() const noexcept {
            double width = Line(top_left, top_right).length();
            double height = Line(top_left, bottom_left).length();
            return width * height;
        }

        double perimeter() const noexcept {
            double w = Line(top_left, top_right).length();
            double h = Line(top_left, bottom_left).length();
            return 2 * (w + h);
        }

        const Point &getTopLeft() const noexcept { return top_left; }
        const Point &getTopRight() const noexcept { return top_right; }
        const Point &getBottomLeft() const noexcept { return bottom_left; }
        const Point &getBottomRight() const noexcept { return bottom_right; }

      private:
        Point top_left;
        Point top_right;
        Point bottom_left;
        Point bottom_right;
    };

} // namespace concord

#pragma once

#include "types_basic.hpp"
#include <cmath>

namespace concord {

    class Circle {
      public:
        Circle() = default;
        Circle(const Point &center_, double radius_) : center(center_), radius(radius_) {}

        double area() const noexcept { return M_PI * radius * radius; }
        double circumference() const noexcept { return 2 * M_PI * radius; }

        const Point &getCenter() const noexcept { return center; }
        void setCenter(const Point &c) noexcept { center = c; }

        bool contains(const Point &p) const noexcept {
            return std::sqrt(std::pow(p.enu.x - center.enu.x, 2) + std::pow(p.enu.y - center.enu.y, 2)) < radius;
        }

        std::vector<Point> as_polygon(int n = 100) const noexcept {
            std::vector<Point> points;
            double theta = 2 * M_PI / n;
            for (int i = 0; i < n; i++) {
                double x = center.enu.x + radius * std::cos(theta * i);
                double y = center.enu.y + radius * std::sin(theta * i);
                Point p;
                p.enu.x = x;
                p.enu.y = y;
                points.push_back(p);
            }
            return points;
        }

        double getRadius() const noexcept { return radius; }
        void setRadius(double r) noexcept { radius = r; }

      private:
        Point center;
        double radius;
    };

} // namespace concord

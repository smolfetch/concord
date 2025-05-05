#pragma once

#include "types_basic.hpp"
#include "types_line.hpp"
#include <cmath>
#include <cstddef>
#include <vector>

namespace concord {

    class Polygon {
      public:
        Polygon() = default;
        explicit Polygon(const std::vector<Point> &pts) : points(pts) {}

        void addPoint(const Point &p) { points.emplace_back(p); }
        void addPoint(const ENU &e, Datum d) { points.emplace_back(Point(e, d)); }
        void addPoint(const WGS &w, Datum d) { points.emplace_back(Point(w, d)); }

        std::size_t numVertices() const noexcept { return points.size(); }
        bool isConnected() const noexcept { return points.size() >= 3; }

        double perimeter() const noexcept {
            if (points.size() < 2)
                return 0.0;
            double per = 0.0;
            for (std::size_t i = 1; i < points.size(); ++i)
                per += Line(points[i - 1], points[i]).length();
            per += Line(points.back(), points.front()).length();
            return per;
        }

        double area() const noexcept {
            if (points.size() < 3)
                return 0.0;
            double a = 0.0;
            for (std::size_t i = 0, j = points.size() - 1; i < points.size(); j = i++) {
                const auto &pi = points[i].enu;
                const auto &pj = points[j].enu;
                a += (pj.x + pi.x) * (pj.y - pi.y);
            }
            return std::abs(a * 0.5);
        }

        bool contains(const Point &p) const noexcept {
            if (points.size() < 3)
                return false;
            bool c = false;
            for (std::size_t i = 0, j = points.size() - 1; i < points.size(); j = i++) {
                const auto &pi = points[i].enu;
                const auto &pj = points[j].enu;
                if (((pi.y > p.enu.y) != (pj.y > p.enu.y)) &&
                    (p.enu.x < (pj.x - pi.x) * (p.enu.y - pi.y) / (pj.y - pi.y) + pi.x))
                    c = !c;
            }
            return c;
        }

        Polygon from_rectangle(const float width, const float height, Datum d = {},
                               Size inflate = Size(1.0, 1.0, 1.0)) const {
            Polygon p;
            p.addPoint(Point(ENU(width * inflate.x / 2.0, height * inflate.y / 2.0, 0.0), d));
            p.addPoint(Point(ENU(width * inflate.x / 2.0, -height * inflate.y / 2.0, 0.0), d));
            p.addPoint(Point(ENU(-width * inflate.x / 2.0, -height * inflate.y / 2.0, 0.0), d));
            p.addPoint(Point(ENU(-width * inflate.x / 2.0, height * inflate.y / 2.0, 0.0), d));
            return p;
        }

        Polygon from_rectangle(Size size, Datum d = {}, Size inflate = Size(1.0, 1.0, 1.0)) const {
            return from_rectangle(size.x, size.y, d, inflate);
        }

        auto begin() noexcept { return points.begin(); }
        auto end() noexcept { return points.end(); }
        auto begin() const noexcept { return points.begin(); }
        auto end() const noexcept { return points.end(); }

      private:
        std::vector<Point> points;
    };

} // namespace concord

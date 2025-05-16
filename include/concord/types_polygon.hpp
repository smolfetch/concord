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

        void from_wgs(std::vector<WGS> pts, Datum d = {}) {
            for (auto &pt : pts) {
                addPoint(Point(pt, d));
            }
        }

        void from_enu(std::vector<ENU> pts, Datum d = {}) {
            for (auto &pt : pts) {
                addPoint(Point(pt, d));
            }
        }

        Polygon from_rectangle(Size size, Datum d = {}, Size inflate = Size(1.0, 1.0, 1.0)) const {
            return from_rectangle(size.x, size.y, d, inflate);
        }

        concord::Bound get_obb(concord::Datum d = {}) const {
            if (points.empty()) {
                return concord::Bound();
            }
            // 1) compute AABB in ENU
            double minX = points[0].enu.x;
            double maxX = minX;
            double minY = points[0].enu.y;
            double maxY = minY;

            for (std::size_t i = 1; i < points.size(); ++i) {
                const auto &e = points[i].enu;
                minX = std::min(minX, e.x);
                maxX = std::max(maxX, e.x);
                minY = std::min(minY, e.y);
                maxY = std::max(maxY, e.y);
            }

            // 2) compute centroid-based “orientation” (same as before)
            double sumX = 0.0, sumY = 0.0;
            for (const auto &p : points) {
                sumX += p.enu.x;
                sumY += p.enu.y;
            }
            double centroidX = sumX / points.size();
            double centroidY = sumY / points.size();
            const auto &first = points[0].enu;
            double orientation_rad = std::atan2(centroidY - first.y, centroidX - first.x);

            // 3) build the concord::Pose
            double centerX = 0.5 * (minX + maxX);
            double centerY = 0.5 * (minY + maxY);
            concord::ENU center_enu(centerX, centerY, 0.0);
            concord::Point center_pt(center_enu, d);

            // Euler takes (roll, pitch, yaw)
            concord::Euler euler(0.0, 0.0, orientation_rad);
            concord::Pose pose(center_pt, euler);

            // 4) build concord::Size (width, height, depth)
            double width = maxX - minX;
            double height = maxY - minY;
            concord::Size size(width, height, 0.0);

            // 5) return the Bound
            return concord::Bound(pose, size);
        }

        auto begin() noexcept { return points.begin(); }
        auto end() noexcept { return points.end(); }
        auto begin() const noexcept { return points.begin(); }
        auto end() const noexcept { return points.end(); }

      private:
        std::vector<Point> points;
    };

} // namespace concord

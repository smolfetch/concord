#pragma once

#include "../core/types_basic.hpp"
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

        Polygon from_vector(std::vector<Point> pts) {
            Polygon p;
            for (auto &pt : pts) {
                p.addPoint(pt);
            }
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

        Bound get_obb(concord::Datum d = {}) const {
            if (points.empty()) {
                return Bound();
            }

            // --- 1) Compute centroid-based orientation (same as before) ---
            double sumX = 0.0, sumY = 0.0;
            for (const auto &p : points) {
                sumX += p.enu.x;
                sumY += p.enu.y;
            }
            double centroidX = sumX / points.size();
            double centroidY = sumY / points.size();

            const auto &first = points[0].enu;
            double orientation_rad = std::atan2(centroidY - first.y, centroidX - first.x);
            double cosO = std::cos(orientation_rad);
            double sinO = std::sin(orientation_rad);

            // --- 2) Rotate all points into this frame and track min/max ---
            double minRotX = std::numeric_limits<double>::infinity();
            double maxRotX = -std::numeric_limits<double>::infinity();
            double minRotY = std::numeric_limits<double>::infinity();
            double maxRotY = -std::numeric_limits<double>::infinity();

            for (const auto &p : points) {
                double x = p.enu.x;
                double y = p.enu.y;
                // rotate (x,y) by +orientation_rad
                double rotX = x * cosO + y * sinO;
                double rotY = -x * sinO + y * cosO;

                minRotX = std::min(minRotX, rotX);
                maxRotX = std::max(maxRotX, rotX);
                minRotY = std::min(minRotY, rotY);
                maxRotY = std::max(maxRotY, rotY);
            }

            // --- 3) Compute width/height in rotated frame ---
            double width = maxRotX - minRotX;
            double height = maxRotY - minRotY;

            // --- 4) Compute center in rotated coords, then un-rotate back ---
            double centerRotX = 0.5 * (minRotX + maxRotX);
            double centerRotY = 0.5 * (minRotY + maxRotY);

            // inverse rotation by -orientation_rad
            double centerX = centerRotX * cosO - centerRotY * sinO;
            double centerY = centerRotX * sinO + centerRotY * cosO;

            // --- 5) Build and return the final Bound ---
            concord::ENU center_enu(centerX, centerY, 0.0);
            concord::Point center_pt(center_enu, d);
            concord::Euler euler(0.0, 0.0, orientation_rad);
            concord::Pose pose(center_pt, euler);

            concord::Size size(width, height, 0.0);
            return concord::Bound(pose, size);
        }

        auto begin() noexcept { return points.begin(); }
        auto end() noexcept { return points.end(); }
        auto begin() const noexcept { return points.begin(); }
        auto end() const noexcept { return points.end(); }

        const std::vector<Point> &getPoints() const noexcept { return points; }

      private:
        std::vector<Point> points;
    };

} // namespace concord

#pragma once

#include "pose.hpp"
#include "size.hpp"
#include "point.hpp"
#include "datum.hpp"
#include <vector>
#include <cmath>

namespace concord {

    struct Bound {
        Pose pose;
        Size size;

        Bound() = default;
        Bound(const Pose &p, const Size &s) : pose(p), size(s) {}
        inline bool is_set() const { return pose.is_set() && size.is_set(); }
        
        // Geometric properties
        inline double volume() const { return size.volume(); }
        inline double area() const { return size.area_xy(); }
        inline Point center() const { return pose.point; }
        
        // Containment testing
        inline bool contains(const Point& point) const {
            Point local_point = pose.inverse_transform_point(point);
            double half_x = size.x * 0.5;
            double half_y = size.y * 0.5;
            double half_z = size.z * 0.5;
            return (std::abs(local_point.enu.x) <= half_x &&
                   std::abs(local_point.enu.y) <= half_y &&
                   std::abs(local_point.enu.z) <= half_z);
        }
        
        // Intersection testing
        inline bool intersects(const Bound& other) const {
            // Simple OBB intersection test using separating axis theorem
            // This is a simplified version - full SAT would check all axes
            std::vector<Point> corners1 = get_corners();
            std::vector<Point> corners2 = other.get_corners();
            
            // Check if any corner of one bound is inside the other
            for (const auto& corner : corners1) {
                if (other.contains(corner)) return true;
            }
            for (const auto& corner : corners2) {
                if (contains(corner)) return true;
            }
            return false;
        }
        
        // Expand bound to include point
        inline void expand_to_include(const Point& point) {
            if (!contains(point)) {
                // Transform point to local coordinates
                Point local_point = pose.inverse_transform_point(point);
                
                // Expand size if necessary
                double half_x = size.x * 0.5;
                double half_y = size.y * 0.5;
                double half_z = size.z * 0.5;
                
                if (std::abs(local_point.enu.x) > half_x) {
                    size.x = std::abs(local_point.enu.x) * 2.0;
                }
                if (std::abs(local_point.enu.y) > half_y) {
                    size.y = std::abs(local_point.enu.y) * 2.0;
                }
                if (std::abs(local_point.enu.z) > half_z) {
                    size.z = std::abs(local_point.enu.z) * 2.0;
                }
            }
        }

        inline std::vector<Point> get_corners(Datum datum = {}) const {
            return pose.get_corners(size, datum);
        }
    };

} // namespace concord

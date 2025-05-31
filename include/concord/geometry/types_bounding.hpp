#pragma once

#include "../core/types_basic.hpp"
#include "../math/types_math.hpp"
#include <array>
#include <limits>
#include <algorithm>

namespace concord {

    // Axis-Aligned Bounding Box (AABB)
    struct AABB {
        Point min_point;
        Point max_point;
        
        AABB() = default;
        AABB(const Point& min_p, const Point& max_p) : min_point(min_p), max_point(max_p) {}
        
        // Create AABB from collection of points
        template<typename Container>
        static AABB fromPoints(const Container& points, const Datum& /* datum */ = {}) {
            if (points.empty()) {
                return AABB{};
            }
            
            auto it = points.begin();
            AABB result;
            result.min_point = *it;
            result.max_point = *it;
            
            for (++it; it != points.end(); ++it) {
                result.expand(*it);
            }
            return result;
        }
        
        inline void expand(const Point& point) {
            min_point.enu.x = std::min(min_point.enu.x, point.enu.x);
            min_point.enu.y = std::min(min_point.enu.y, point.enu.y);
            min_point.enu.z = std::min(min_point.enu.z, point.enu.z);
            
            max_point.enu.x = std::max(max_point.enu.x, point.enu.x);
            max_point.enu.y = std::max(max_point.enu.y, point.enu.y);
            max_point.enu.z = std::max(max_point.enu.z, point.enu.z);
        }
        
        inline bool contains(const Point& point) const {
            return (point.enu.x >= min_point.enu.x && point.enu.x <= max_point.enu.x &&
                    point.enu.y >= min_point.enu.y && point.enu.y <= max_point.enu.y &&
                    point.enu.z >= min_point.enu.z && point.enu.z <= max_point.enu.z);
        }
        
        inline bool intersects(const AABB& other) const {
            return (min_point.enu.x <= other.max_point.enu.x && max_point.enu.x >= other.min_point.enu.x &&
                    min_point.enu.y <= other.max_point.enu.y && max_point.enu.y >= other.min_point.enu.y &&
                    min_point.enu.z <= other.max_point.enu.z && max_point.enu.z >= other.min_point.enu.z);
        }
        
        inline Point center() const {
            Point center_pt;
            center_pt.enu.x = (min_point.enu.x + max_point.enu.x) * 0.5;
            center_pt.enu.y = (min_point.enu.y + max_point.enu.y) * 0.5;
            center_pt.enu.z = (min_point.enu.z + max_point.enu.z) * 0.5;
            return center_pt;
        }
        
        inline Size size() const {
            return Size{
                max_point.enu.x - min_point.enu.x,
                max_point.enu.y - min_point.enu.y,
                max_point.enu.z - min_point.enu.z
            };
        }
        
        inline double volume() const {
            auto s = size();
            return s.x * s.y * s.z;
        }
        
        inline double surface_area() const {
            auto s = size();
            return 2.0 * (s.x * s.y + s.y * s.z + s.z * s.x);
        }
        
        inline std::array<Point, 8> corners() const {
            return {
                Point(ENU{min_point.enu.x, min_point.enu.y, min_point.enu.z}, Datum{}),
                Point(ENU{max_point.enu.x, min_point.enu.y, min_point.enu.z}, Datum{}),
                Point(ENU{max_point.enu.x, max_point.enu.y, min_point.enu.z}, Datum{}),
                Point(ENU{min_point.enu.x, max_point.enu.y, min_point.enu.z}, Datum{}),
                Point(ENU{min_point.enu.x, min_point.enu.y, max_point.enu.z}, Datum{}),
                Point(ENU{max_point.enu.x, min_point.enu.y, max_point.enu.z}, Datum{}),
                Point(ENU{max_point.enu.x, max_point.enu.y, max_point.enu.z}, Datum{}),
                Point(ENU{min_point.enu.x, max_point.enu.y, max_point.enu.z}, Datum{})
            };
        }
    };
    
    // Oriented Bounding Box (OBB) - more precise than AABB
    struct OBB {
        Point center;
        Size half_extents;
        Euler orientation;
        
        OBB() = default;
        OBB(const Point& c, const Size& he, const Euler& orient) 
            : center(c), half_extents(he), orientation(orient) {}
        
        // Create OBB from points using PCA or convex hull + rotating calipers
        template<typename Container>
        static OBB fromPoints(const Container& points, const Datum& datum = {}) {
            if (points.empty()) {
                return OBB{};
            }
            
            // For now, use AABB as approximation
            // TODO: Implement proper PCA-based OBB fitting
            auto aabb = AABB::fromPoints(points, datum);
            return OBB{
                aabb.center(),
                Size{aabb.size().x * 0.5, aabb.size().y * 0.5, aabb.size().z * 0.5},
                Euler{}
            };
        }
        
        inline bool contains(const Point& point) const {
            // Transform point to OBB local space
            double dx = point.enu.x - center.enu.x;
            double dy = point.enu.y - center.enu.y;
            double dz = point.enu.z - center.enu.z;
            
            // Apply inverse rotation
            double cos_yaw = std::cos(-orientation.yaw);
            double sin_yaw = std::sin(-orientation.yaw);
            // double cos_pitch = std::cos(-orientation.pitch);  // Currently unused
            // double sin_pitch = std::sin(-orientation.pitch);  // Currently unused
            // double cos_roll = std::cos(-orientation.roll);    // Currently unused
            // double sin_roll = std::sin(-orientation.roll);    // Currently unused
            
            // Simplified rotation (yaw only for 2D case)
            double local_x = dx * cos_yaw - dy * sin_yaw;
            double local_y = dx * sin_yaw + dy * cos_yaw;
            double local_z = dz;
            
            return (std::abs(local_x) <= half_extents.x &&
                    std::abs(local_y) <= half_extents.y &&
                    std::abs(local_z) <= half_extents.z);
        }
        
        inline std::array<Point, 8> corners(const Datum& datum = {}) const {
            std::array<Point, 8> corners;
            
            // Local corners relative to center
            std::array<std::array<double, 3>, 8> local_corners = {{
                {{-half_extents.x, -half_extents.y, -half_extents.z}},
                {{+half_extents.x, -half_extents.y, -half_extents.z}},
                {{+half_extents.x, +half_extents.y, -half_extents.z}},
                {{-half_extents.x, +half_extents.y, -half_extents.z}},
                {{-half_extents.x, -half_extents.y, +half_extents.z}},
                {{+half_extents.x, -half_extents.y, +half_extents.z}},
                {{+half_extents.x, +half_extents.y, +half_extents.z}},
                {{-half_extents.x, +half_extents.y, +half_extents.z}}
            }};
            
            double cos_yaw = std::cos(orientation.yaw);
            double sin_yaw = std::sin(orientation.yaw);
            
            for (size_t i = 0; i < 8; ++i) {
                // Apply rotation and translation
                double lx = local_corners[i][0];
                double ly = local_corners[i][1];
                double lz = local_corners[i][2];
                
                double world_x = center.enu.x + (lx * cos_yaw - ly * sin_yaw);
                double world_y = center.enu.y + (lx * sin_yaw + ly * cos_yaw);
                double world_z = center.enu.z + lz;
                
                corners[i] = Point{ENU{world_x, world_y, world_z}, datum};
            }
            
            return corners;
        }
    };
    
    // Sphere/Circle bounding volume
    struct BoundingSphere {
        Point center;
        double radius = 0.0;
        
        BoundingSphere() = default;
        BoundingSphere(const Point& c, double r) : center(c), radius(r) {}
        
        template<typename Container>
        static BoundingSphere fromPoints(const Container& points) {
            if (points.empty()) {
                return BoundingSphere{};
            }
            
            // Naive algorithm: use centroid and max distance
            Point centroid;
            double sum_x = 0, sum_y = 0, sum_z = 0;
            size_t count = 0;
            
            for (const auto& point : points) {
                sum_x += point.enu.x;
                sum_y += point.enu.y;
                sum_z += point.enu.z;
                ++count;
            }
            
            centroid.enu.x = sum_x / count;
            centroid.enu.y = sum_y / count;
            centroid.enu.z = sum_z / count;
            
            double max_dist = 0.0;
            for (const auto& point : points) {
                double dx = point.enu.x - centroid.enu.x;
                double dy = point.enu.y - centroid.enu.y;
                double dz = point.enu.z - centroid.enu.z;
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                max_dist = std::max(max_dist, dist);
            }
            
            return BoundingSphere{centroid, max_dist};
        }
        
        inline bool contains(const Point& point) const {
            double dx = point.enu.x - center.enu.x;
            double dy = point.enu.y - center.enu.y;
            double dz = point.enu.z - center.enu.z;
            double dist_sq = dx*dx + dy*dy + dz*dz;
            return dist_sq <= radius * radius;
        }
        
        inline bool intersects(const BoundingSphere& other) const {
            double dx = other.center.enu.x - center.enu.x;
            double dy = other.center.enu.y - center.enu.y;
            double dz = other.center.enu.z - center.enu.z;
            double dist_sq = dx*dx + dy*dy + dz*dz;
            double sum_radii = radius + other.radius;
            return dist_sq <= sum_radii * sum_radii;
        }
        
        inline double volume() const {
            return (4.0/3.0) * M_PI * radius * radius * radius;
        }
        
        inline double surface_area() const {
            return 4.0 * M_PI * radius * radius;
        }
    };

} // namespace concord

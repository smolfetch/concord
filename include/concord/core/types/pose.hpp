#pragma once

#include "point.hpp"
#include "euler.hpp"
#include "quaternion.hpp"
#include "size.hpp"
#include "enu.hpp"
#include "wgs.hpp"
#include "datum.hpp"
#include "../../math/types_math.hpp"
#include <vector>
#include <array>
#include <cmath>

namespace concord {

    struct Pose {
        Point point;
        Euler angle;

        Pose() = default;
        Pose(const Point &p, const Euler &a) : point(p), angle(a) {}
        Pose(float x, float y, float yaw)
            : point(Point{ENU{x, y, 0.0f}, WGS{0.0f, 0.0f, 0.0f}}), angle(Euler{0.0f, 0.0f, yaw}) {}
        explicit Pose(const Point &p, const Quaternion &q) noexcept : point(p), angle(q) {}
        inline bool is_set() const { return point.is_set() && angle.is_set(); }
        
        // Transformation operations
        inline Point transform_point(const Point& local_point) const {
            Vec3d local_vec = local_point.enu.to_vec3();
            Vec3d rotated = angle.rotate(local_vec);
            ENU transformed_enu = point.enu + ENU::from_vec3(rotated);
            return Point{transformed_enu, transformed_enu.toWGS(Datum{})};
        }
        
        inline Point inverse_transform_point(const Point& world_point) const {
            ENU relative_enu = world_point.enu - point.enu;
            Vec3d relative_vec = relative_enu.to_vec3();
            Euler inverse_angle = angle * -1.0;
            Vec3d local_vec = inverse_angle.rotate(relative_vec);
            ENU local_enu = ENU::from_vec3(local_vec);
            return Point{local_enu, local_enu.toWGS(Datum{})};
        }
        
        // Pose composition
        inline Pose operator*(const Pose& other) const {
            Point transformed_point = transform_point(other.point);
            Euler combined_angle = angle + other.angle;
            return Pose{transformed_point, combined_angle.normalized()};
        }
        
        // Inverse pose
        inline Pose inverse() const {
            Euler inv_angle = angle * -1.0;
            Vec3d neg_pos = point.enu.to_vec3() * -1.0;
            Vec3d rotated_neg_pos = inv_angle.rotate(neg_pos);
            ENU inv_enu = ENU::from_vec3(rotated_neg_pos);
            Point inv_point{inv_enu, inv_enu.toWGS(Datum{})};
            return Pose{inv_point, inv_angle};
        }

        inline std::vector<Point> get_corners(Size size, Datum datum = {}) const {
            std::vector<Point> points;
            // precompute
            double c = std::cos(angle.yaw);
            double s = std::sin(angle.yaw);
            double halfW = size.x * 0.5;
            double halfH = size.y * 0.5;
            double cx = point.enu.x;
            double cy = point.enu.y;
            // local corners (in CCW order, for example)
            std::array<std::pair<double, double>, 4> local = {
                {{+halfW, +halfH}, {+halfW, -halfH}, {-halfW, -halfH}, {-halfW, +halfH}}};
            for (auto [lx, ly] : local) {
                // rotate
                double rx = c * lx - s * ly;
                double ry = s * lx + c * ly;
                // translate back into world‐ENU
                double wx = cx + rx;
                double wy = cy + ry;
                points.emplace_back(Point(ENU{wx, wy, 0.0}, datum));
            }
            return points;
        }
    };

} // namespace concord

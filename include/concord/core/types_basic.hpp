#pragma once

#include "wgs_to_enu.hpp"
#include "../math/types_math.hpp"
#include "../errors/error_handling.hpp"
#include <array>
#include <tuple>
#include <variant>
#include <vector>
#include <cmath>

namespace concord {

    enum class CRS { ENU, WGS };

    struct ENU;        // forward declaration
    struct WGS;        // forward declaration
    struct Euler;      // forward declaration
    struct Quaternion; // forward declaration

    struct Datum {
        double lat = 0.0;
        double lon = 0.0;
        double alt = 0.0;

        inline bool is_set() const { return lat != 0.0 && lon != 0.0; }
    };

    struct WGS {
        double lat = 0.0;
        double lon = 0.0;
        double alt = 0.0;

        WGS(double lat_, double lon_, double alt_) : lat(lat_), lon(lon_), alt(alt_) {
            validation::validate_latitude(lat);
            validation::validate_longitude(lon);
            validation::validate_altitude(alt);
        }

        WGS() = default;
        inline ENU toENU(const Datum &datum) const;
        inline operator Datum() const noexcept { return Datum{lat, lon, alt}; }
        inline bool is_set() const { return lat != 0.0 && lon != 0.0; }
        
        // Mathematical operations (for small displacements)
        inline WGS operator+(const WGS& offset) const { 
            return WGS{lat + offset.lat, lon + offset.lon, alt + offset.alt}; 
        }
        inline WGS operator-(const WGS& offset) const { 
            return WGS{lat - offset.lat, lon - offset.lon, alt - offset.alt}; 
        }
        
        // Distance calculation using haversine formula
        inline double distance_to(const WGS& other) const {
            validation::validate_finite(lat, "latitude");
            validation::validate_finite(lon, "longitude");
            validation::validate_finite(other.lat, "other latitude");
            validation::validate_finite(other.lon, "other longitude");
            
            const double R = 6371000.0; // Earth radius in meters
            double dlat = (other.lat - lat) * M_PI / 180.0;
            double dlon = (other.lon - lon) * M_PI / 180.0;
            double a = std::sin(dlat/2) * std::sin(dlat/2) +
                      std::cos(lat * M_PI / 180.0) * std::cos(other.lat * M_PI / 180.0) *
                      std::sin(dlon/2) * std::sin(dlon/2);
            double c = 2 * safe_math::safe_asin(safe_math::safe_sqrt(a, "haversine"), "haversine");
            return R * c;
        }
        
        // Bearing calculation
        inline double bearing_to(const WGS& other) const {
            double dlon = (other.lon - lon) * M_PI / 180.0;
            double lat1_rad = lat * M_PI / 180.0;
            double lat2_rad = other.lat * M_PI / 180.0;
            
            double y = std::sin(dlon) * std::cos(lat2_rad);
            double x = std::cos(lat1_rad) * std::sin(lat2_rad) - 
                      std::sin(lat1_rad) * std::cos(lat2_rad) * std::cos(dlon);
            double bearing = std::atan2(y, x) * 180.0 / M_PI;
            return std::fmod(bearing + 360.0, 360.0);
        }
    };

    struct ENU {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        ENU(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

        ENU() = default;
        inline WGS toWGS(const Datum &datum) const;
        inline bool is_set() const { return x != 0.0 && y != 0.0; }
        
        // Mathematical operations
        inline ENU operator+(const ENU& other) const { return ENU{x + other.x, y + other.y, z + other.z}; }
        inline ENU operator-(const ENU& other) const { return ENU{x - other.x, y - other.y, z - other.z}; }
        inline ENU operator*(double scale) const { return ENU{x * scale, y * scale, z * scale}; }
        inline ENU operator/(double scale) const { return ENU{x / scale, y / scale, z / scale}; }
        
        inline ENU& operator+=(const ENU& other) { x += other.x; y += other.y; z += other.z; return *this; }
        inline ENU& operator-=(const ENU& other) { x -= other.x; y -= other.y; z -= other.z; return *this; }
        inline ENU& operator*=(double scale) { x *= scale; y *= scale; z *= scale; return *this; }
        inline ENU& operator/=(double scale) { x /= scale; y /= scale; z /= scale; return *this; }
        
        // Distance and magnitude operations
        inline double magnitude() const { return std::sqrt(x * x + y * y + z * z); }
        inline double distance_to(const ENU& other) const { 
            return std::sqrt((x - other.x) * (x - other.x) + 
                           (y - other.y) * (y - other.y) + 
                           (z - other.z) * (z - other.z)); 
        }
        inline double distance_to_2d(const ENU& other) const {
            return std::sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
        }
        
        // Conversion to Vec3d
        inline Vec3d to_vec3() const { return Vec3d{x, y, z}; }
        inline static ENU from_vec3(const Vec3d& v) { return ENU{v[0], v[1], v[2]}; }
    };

    inline ENU WGS::toENU(const Datum &datum) const {
        auto enu = gps_to_enu(lat, lon, alt, datum.lat, datum.lon, datum.alt);
        return ENU{std::get<0>(enu), std::get<1>(enu), std::get<2>(enu)};
    }

    inline WGS ENU::toWGS(const Datum &datum) const {
        auto wgs = enu_to_gps(x, y, z, datum.lat, datum.lon, datum.alt);
        return WGS{std::get<0>(wgs), std::get<1>(wgs), std::get<2>(wgs)};
    }

    struct Quaternion {
        double w = 0.0;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        Quaternion() = default;
        Quaternion(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}
        explicit Quaternion(const Euler &e) noexcept;
        inline bool is_set() const { return w != 0.0 || x != 0.0 || y != 0.0 || z != 0.0; }
        
        // Mathematical operations
        inline Quaternion operator+(const Quaternion& other) const {
            return Quaternion{w + other.w, x + other.x, y + other.y, z + other.z};
        }
        inline Quaternion operator-(const Quaternion& other) const {
            return Quaternion{w - other.w, x - other.x, y - other.y, z - other.z};
        }
        inline Quaternion operator*(const Quaternion& other) const {
            return Quaternion{
                w * other.w - x * other.x - y * other.y - z * other.z,
                w * other.x + x * other.w + y * other.z - z * other.y,
                w * other.y - x * other.z + y * other.w + z * other.x,
                w * other.z + x * other.y - y * other.x + z * other.w
            };
        }
        inline Quaternion operator*(double scale) const {
            return Quaternion{w * scale, x * scale, y * scale, z * scale};
        }
        
        // Quaternion operations
        inline double norm() const { 
            validation::validate_finite(w, "quaternion w");
            validation::validate_finite(x, "quaternion x");
            validation::validate_finite(y, "quaternion y");
            validation::validate_finite(z, "quaternion z");
            return std::sqrt(w*w + x*x + y*y + z*z); 
        }
        
        inline Quaternion normalized() const {
            double n = norm();
            if (n < 1e-15) {
                throw MathematicalException("cannot normalize zero quaternion");
            }
            return (*this) * (1.0 / n);
        }
        
        inline Quaternion conjugate() const { return Quaternion{w, -x, -y, -z}; }
        
        inline Quaternion inverse() const {
            double n2 = w*w + x*x + y*y + z*z;
            if (n2 < 1e-15) {
                throw MathematicalException("cannot invert zero quaternion");
            }
            return conjugate() * (1.0 / n2);
        }
        
        // Rotation operations
        inline Vec3d rotate(const Vec3d& v) const {
            Quaternion qv{0, v[0], v[1], v[2]};
            Quaternion result = (*this) * qv * conjugate();
            return Vec3d{result.x, result.y, result.z};
        }
        
        // SLERP interpolation
        inline static Quaternion slerp(const Quaternion& q1, const Quaternion& q2, double t) {
            double dot = q1.w * q2.w + q1.x * q2.x + q1.y * q2.y + q1.z * q2.z;
            if (dot < 0.0) {
                return slerp(q1, Quaternion{-q2.w, -q2.x, -q2.y, -q2.z}, t);
            }
            if (dot > 0.9995) {
                // Linear interpolation for very close quaternions
                Quaternion result = q1 + (q2 - q1) * t;
                return result.normalized();
            }
            double theta = std::acos(std::abs(dot));
            double sin_theta = std::sin(theta);
            double t1 = std::sin((1.0 - t) * theta) / sin_theta;
            double t2 = std::sin(t * theta) / sin_theta;
            return q1 * t1 + q2 * t2;
        }
    };

    struct Euler {
        double roll = 0.0;  // rotation about x-axis
        double pitch = 0.0; // rotation about y-axis
        double yaw = 0.0;   // rotation about z-axis

        Euler() = default;
        Euler(double roll_, double pitch_, double yaw_) : roll(roll_), pitch(pitch_), yaw(yaw_) {}
        explicit Euler(const Quaternion &q) noexcept;
        inline bool is_set() const { return roll != 0.0 || pitch != 0.0 || yaw != 0.0; }

        inline double yaw_cos() const { return std::cos(yaw * 0.5); }
        inline double yaw_sin() const { return std::sin(yaw * 0.5); }
        
        // Mathematical operations
        inline Euler operator+(const Euler& other) const {
            return Euler{roll + other.roll, pitch + other.pitch, yaw + other.yaw};
        }
        inline Euler operator-(const Euler& other) const {
            return Euler{roll - other.roll, pitch - other.pitch, yaw - other.yaw};
        }
        inline Euler operator*(double scale) const {
            return Euler{roll * scale, pitch * scale, yaw * scale};
        }
        
        // Angle normalization
        inline Euler normalized() const {
            auto normalize_angle = [](double angle) {
                while (angle > M_PI) angle -= 2.0 * M_PI;
                while (angle < -M_PI) angle += 2.0 * M_PI;
                return angle;
            };
            return Euler{normalize_angle(roll), normalize_angle(pitch), normalize_angle(yaw)};
        }
        
        // Convert to rotation matrix (3x3)
        inline Mat3d to_rotation_matrix() const {
            double cr = std::cos(roll), sr = std::sin(roll);
            double cp = std::cos(pitch), sp = std::sin(pitch);
            double cy = std::cos(yaw), sy = std::sin(yaw);
            
            Mat3d R;
            R[0][0] = cy * cp;
            R[0][1] = cy * sp * sr - sy * cr;
            R[0][2] = cy * sp * cr + sy * sr;
            R[1][0] = sy * cp;
            R[1][1] = sy * sp * sr + cy * cr;
            R[1][2] = sy * sp * cr - cy * sr;
            R[2][0] = -sp;
            R[2][1] = cp * sr;
            R[2][2] = cp * cr;
            return R;
        }
        
        // Rotate a vector
        inline Vec3d rotate(const Vec3d& v) const {
            Mat3d R = to_rotation_matrix();
            return R * v;
        }
    };

    // — Definitions —

    inline Euler::Euler(const Quaternion &q) noexcept {
        // roll (x-axis rotation)
        double sinr_cosp = 2.0 * (q.w * q.x + q.y * q.z);
        double cosr_cosp = 1.0 - 2.0 * (q.x * q.x + q.y * q.y);
        roll = std::atan2(sinr_cosp, cosr_cosp);

        // pitch (y-axis rotation)
        double sinp = 2.0 * (q.w * q.y - q.z * q.x);
        if (std::abs(sinp) >= 1.0) {
            pitch = std::copysign(M_PI / 2.0, sinp);
        } else {
            pitch = std::asin(sinp);
        }

        // yaw (z-axis rotation)
        double siny_cosp = 2.0 * (q.w * q.z + q.x * q.y);
        double cosy_cosp = 1.0 - 2.0 * (q.y * q.y + q.z * q.z);
        yaw = std::atan2(siny_cosp, cosy_cosp);
    }

    inline Quaternion::Quaternion(const Euler &e) noexcept {
        double cy = std::cos(e.yaw * 0.5);
        double sy = std::sin(e.yaw * 0.5);
        double cp = std::cos(e.pitch * 0.5);
        double sp = std::sin(e.pitch * 0.5);
        double cr = std::cos(e.roll * 0.5);
        double sr = std::sin(e.roll * 0.5);

        w = cr * cp * cy + sr * sp * sy;
        x = sr * cp * cy - cr * sp * sy;
        y = cr * sp * cy + sr * cp * sy;
        z = cr * cp * sy - sr * sp * cy;
    }

    struct Size {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        Size() = default;
        Size(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
        Size(double s) : x(s), y(s), z(s) {} // Uniform scaling
        inline bool is_set() const { return x != 0.0 || y != 0.0 || z != 0.0; }
        
        // Mathematical operations
        inline Size operator+(const Size& other) const { return Size{x + other.x, y + other.y, z + other.z}; }
        inline Size operator-(const Size& other) const { return Size{x - other.x, y - other.y, z - other.z}; }
        inline Size operator*(double scale) const { return Size{x * scale, y * scale, z * scale}; }
        inline Size operator/(double scale) const { return Size{x / scale, y / scale, z / scale}; }
        inline Size operator*(const Size& other) const { return Size{x * other.x, y * other.y, z * other.z}; }
        
        // Volume and area calculations
        inline double volume() const { return x * y * z; }
        inline double area_xy() const { return x * y; }
        inline double area_xz() const { return x * z; }
        inline double area_yz() const { return y * z; }
        inline double diagonal() const { return std::sqrt(x*x + y*y + z*z); }
        inline double diagonal_2d() const { return std::sqrt(x*x + y*y); }
        
        // Utility functions
        inline Size abs() const { return Size{std::abs(x), std::abs(y), std::abs(z)}; }
        inline Size max(const Size& other) const { 
            return Size{std::max(x, other.x), std::max(y, other.y), std::max(z, other.z)}; 
        }
        inline Size min(const Size& other) const { 
            return Size{std::min(x, other.x), std::min(y, other.y), std::min(z, other.z)}; 
        }
    };

    struct Point {
        ENU enu;
        WGS wgs;

        Point() = default;
        Point(const ENU &e, const WGS &w) : enu(e), wgs(w) {}
        Point(const ENU &e, Datum d) : enu(e), wgs(e.toWGS(d)) {}
        Point(const WGS &w, Datum d) : enu(w.toENU(d)), wgs(w) {}

        Point(double x, double y, double z = 0.0, const Datum &datum = Datum()) 
            : enu(x, y, z), wgs(enu.toWGS(datum)) {}
            
        inline bool is_set() const { return enu.is_set() && wgs.is_set(); }
        
        // Mathematical operations
        inline Point operator+(const ENU& offset) const { 
            return Point{enu + offset, wgs}; 
        }
        inline Point operator-(const ENU& offset) const { 
            return Point{enu - offset, wgs}; 
        }
        inline double distance_to(const Point& other) const {
            return enu.distance_to(other.enu);
        }
        inline double distance_to_2d(const Point& other) const {
            return enu.distance_to_2d(other.enu);
        }
    };

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

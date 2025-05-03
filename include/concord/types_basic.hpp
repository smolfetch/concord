#pragma once

#include "wgs_to_enu.hpp"
#include <tuple>
#include <variant>
#include <vector>

namespace concord {
    struct ENU;        // forward declaration
    struct WGS;        // forward declaration
    struct Euler;      // forward declaration
    struct Quaternion; // forward declaration

    struct Datum {
        double lat = 0.0;
        double lon = 0.0;
        double alt = 0.0;
    };

    struct WGS {
        double lat = 0.0;
        double lon = 0.0;
        double alt = 0.0;

        WGS(double lat_, double lon_, double alt_) : lat(lat_), lon(lon_), alt(alt_) {}

        WGS() = default;
        inline ENU toENU(const Datum &datum) const;
        operator Datum() const noexcept { return Datum{lat, lon, alt}; }
    };

    struct ENU {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        ENU(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

        ENU() = default;
        inline WGS toWGS(const Datum &datum) const;
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
    };

    struct Euler {
        double roll = 0.0;  // rotation about x-axis
        double pitch = 0.0; // rotation about y-axis
        double yaw = 0.0;   // rotation about z-axis

        Euler() = default;
        Euler(double roll_, double pitch_, double yaw_) : roll(roll_), pitch(pitch_), yaw(yaw_) {}
        explicit Euler(const Quaternion &q) noexcept;
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
    };

    struct Point {
        ENU enu;
        WGS wgs;

        Point() = default;
        Point(const ENU &e, const WGS &w) : enu(e), wgs(w) {}
        Point(const ENU &e, Datum d) : enu(e), wgs(e.toWGS(d)) {}
        Point(const WGS &w, Datum d) : wgs(w), enu(w.toENU(d)) {}

        WGS fromENU(const Datum &d) const { return enu.toWGS(d); }
        ENU toENU(const Datum &d) const { return wgs.toENU(d); }
    };

    struct Pose {
        Point point;
        Euler angle;

        Pose() = default;
        Pose(const Point &p, const Euler &a) : point(p), angle(a) {}
        explicit Pose(const Point &p, const Quaternion &q) noexcept : point(p), angle(q) {}
    };

} // namespace concord

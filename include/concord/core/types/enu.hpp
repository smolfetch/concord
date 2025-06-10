#pragma once

#include "datum.hpp"
#include "wgs.hpp"
#include "../wgs_to_enu.hpp"
#include "../../math/math.hpp"
#include <cmath>
#include <tuple>

namespace concord {

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

    // Implementation of cross-dependencies
    inline ENU WGS::toENU(const Datum &datum) const {
        auto enu = gps_to_enu(lat, lon, alt, datum.lat, datum.lon, datum.alt);
        return ENU{std::get<0>(enu), std::get<1>(enu), std::get<2>(enu)};
    }

    inline WGS ENU::toWGS(const Datum &datum) const {
        auto wgs = enu_to_gps(x, y, z, datum.lat, datum.lon, datum.alt);
        return WGS{std::get<0>(wgs), std::get<1>(wgs), std::get<2>(wgs)};
    }

} // namespace concord

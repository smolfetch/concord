#pragma once

#include "wgs.hpp"
#include "../wgs_to_enu.hpp"
#include <tuple>

namespace concord {

    // ECEF (Earth-Centered, Earth-Fixed) coordinate type
    struct ECEF {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        
        ECEF() = default;
        ECEF(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
        
        inline WGS toWGS() const {
            auto [lat, lon, alt] = ecef_to_gps(x, y, z);
            return WGS{lat, lon, alt};
        }
        
        inline static ECEF fromWGS(const WGS& wgs) {
            auto [x, y, z] = gps_to_ecef(wgs.lat, wgs.lon, wgs.alt);
            return ECEF{x, y, z};
        }
        
        inline bool is_set() const { return x != 0.0 || y != 0.0 || z != 0.0; }
    };

} // namespace concord

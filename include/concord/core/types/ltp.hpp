#pragma once

#include "enu.hpp"

namespace concord {

    // Local Tangent Plane (LTP) coordinate system
    struct LTP {
        double north = 0.0;  // Same as ENU.y
        double east = 0.0;   // Same as ENU.x
        double up = 0.0;     // Same as ENU.z
        
        LTP() = default;
        LTP(double n, double e, double u) : north(n), east(e), up(u) {}
        
        inline ENU toENU() const {
            return ENU{east, north, up};
        }
        
        inline static LTP fromENU(const ENU& enu) {
            return LTP{enu.y, enu.x, enu.z};
        }
        
        inline bool is_set() const { return north != 0.0 || east != 0.0; }
    };

} // namespace concord

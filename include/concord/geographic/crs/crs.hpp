#pragma once

// Individual coordinate system types
#include "../wgs_to_enu.hpp"
#include "../wgs_to_utm.hpp"
#include "ecef.hpp"
#include "enu.hpp"
#include "ltp.hpp"
#include "utm.hpp"
#include "wgs.hpp"

namespace concord {

    // Implementation of WGS::toENU method
    inline ENU WGS::toENU(const Datum &datum) const {
        auto enu = gps_to_enu(lat, lon, alt, datum.lat, datum.lon, datum.alt);
        return ENU{std::get<0>(enu), std::get<1>(enu), std::get<2>(enu), datum};
    }

    // Implementation of ENU::toWGS method
    inline WGS ENU::toWGS() const {
        auto wgs = enu_to_gps(x, y, z, datum.lat, datum.lon, datum.alt);
        return WGS{std::get<0>(wgs), std::get<1>(wgs), std::get<2>(wgs)};
    }

} // namespace concord

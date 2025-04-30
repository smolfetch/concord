#pragma once

#include "wgs_to_enu.hpp"

namespace concord {
    struct ENU;
    struct WGS;

    struct WGS {
        double lat;
        double lon;
        double alt;

        inline ENU toENU(WGS current, WGS datum);
    };

    struct ENU {
        double x;
        double y;
        double z;

        inline WGS toWGS(ENU current, WGS datum);
    };

    inline ENU WGS::toENU(WGS current, WGS datum) {
        std::tuple<double, double, double> enu =
            gps_to_enu(current.lat, current.lon, current.alt, datum.lat, datum.lon, datum.alt);
        return ENU{std::get<0>(enu), std::get<1>(enu), std::get<2>(enu)};
    }

    inline WGS ENU::toWGS(ENU current, WGS datum) {
        std::tuple<double, double, double> wgs =
            enu_to_gps(current.x, current.y, current.z, datum.lat, datum.lon, datum.alt);
        return WGS{std::get<0>(wgs), std::get<1>(wgs), std::get<2>(wgs)};
    }

} // namespace concord

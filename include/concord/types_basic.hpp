#pragma once

#include "wgs_to_enu.hpp"
#include <variant>
#include <vector>

namespace concord {
    struct ENU;
    struct WGS;

    struct Datum {
        double lat;
        double lon;
        double alt;
    };

    struct WGS {
        double lat;
        double lon;
        double alt;

        inline ENU toENU(Datum datum);
        operator Datum() const noexcept { return Datum{lat, lon, alt}; }
    };

    struct ENU {
        double x;
        double y;
        double z;

        inline WGS toWGS(Datum datum);
    };

    inline ENU WGS::toENU(Datum datum) {
        std::tuple<double, double, double> enu = gps_to_enu(lat, lon, alt, datum.lat, datum.lon, datum.alt);
        return ENU{std::get<0>(enu), std::get<1>(enu), std::get<2>(enu)};
    }

    inline WGS ENU::toWGS(Datum datum) {
        std::tuple<double, double, double> wgs = enu_to_gps(x, y, z, datum.lat, datum.lon, datum.alt);
        return WGS{std::get<0>(wgs), std::get<1>(wgs), std::get<2>(wgs)};
    }

    // using Point = std::variant<ENU, WGS>;
    struct Point {
        ENU enu;
        WGS wgs;
    };

} // namespace concord

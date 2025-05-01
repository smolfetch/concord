#pragma once

#include "wgs_to_enu.hpp"
#include <tuple>
#include <variant>
#include <vector>

namespace concord {
    struct ENU; // forward declaration
    struct WGS; // forward declaration

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

        inline ENU toENU(const Datum &datum) const;
        operator Datum() const noexcept { return Datum{lat, lon, alt}; }
    };

    struct ENU {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        ENU(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

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

        Point(const ENU &e, const WGS &w) : enu(e), wgs(w) {}
        explicit Point(const ENU &e, Datum d = {}) : enu(e), wgs(e.toWGS(d)) {}
        explicit Point(const WGS &w, Datum d = {}) : wgs(w), enu(w.toENU(d)) {}
    };

} // namespace concord

#pragma once

#include "wgs_to_enu.hpp"
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

    struct Point {
        ENU enu;
        WGS wgs;
    };

    struct Line {
        Point start;
        Point end;
    };

    struct Path {
        std::vector<Point> points;
    };

    struct Polygon {
        std::vector<Point> points;
        inline bool is_connected() {
            if (points.size() < 3 || points.begin() == points.end()) {
                return false;
            }
            return true;
        }
    };

    struct Circle {
        Point center;
        double radius;
    };

    struct Square {
        Point center;
        double side;
    };

    struct Rectangle {
        Point top_left;
        Point top_right;
        Point bottom_left;
        Point bottom_right;
    };

} // namespace concord

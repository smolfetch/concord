#pragma once

#include "crs/crs.hpp"
#include <cmath>
#include <optional>
#include <string>
#include <tuple>

namespace concord {

    // Coordinate system conversion utilities
    namespace coords {

        // Distance calculations accounting for Earth's curvature
        inline double haversineDistance(const WGS &a, const WGS &b) {
            constexpr double R = 6371000.0; // Earth's radius in meters

            double lat1_rad = a.lat * M_PI / 180.0;
            double lat2_rad = b.lat * M_PI / 180.0;
            double dlat_rad = (b.lat - a.lat) * M_PI / 180.0;
            double dlon_rad = (b.lon - a.lon) * M_PI / 180.0;

            double a_val = std::sin(dlat_rad / 2) * std::sin(dlat_rad / 2) +
                           std::cos(lat1_rad) * std::cos(lat2_rad) * std::sin(dlon_rad / 2) * std::sin(dlon_rad / 2);

            double c = 2 * std::atan2(std::sqrt(a_val), std::sqrt(1 - a_val));

            return R * c;
        }

        // Bearing calculation
        inline double bearing(const WGS &from, const WGS &to) {
            double lat1_rad = from.lat * M_PI / 180.0;
            double lat2_rad = to.lat * M_PI / 180.0;
            double dlon_rad = (to.lon - from.lon) * M_PI / 180.0;

            double y = std::sin(dlon_rad) * std::cos(lat2_rad);
            double x =
                std::cos(lat1_rad) * std::sin(lat2_rad) - std::sin(lat1_rad) * std::cos(lat2_rad) * std::cos(dlon_rad);

            double bearing_rad = std::atan2(y, x);
            return std::fmod(bearing_rad * 180.0 / M_PI + 360.0, 360.0);
        }

        // Point at distance and bearing
        inline WGS pointAtDistanceAndBearing(const WGS &origin, double distance, double bearing_deg) {
            constexpr double R = 6371000.0; // Earth's radius in meters

            double lat1_rad = origin.lat * M_PI / 180.0;
            double lon1_rad = origin.lon * M_PI / 180.0;
            double bearing_rad = bearing_deg * M_PI / 180.0;
            double angular_distance = distance / R;

            double lat2_rad = std::asin(std::sin(lat1_rad) * std::cos(angular_distance) +
                                        std::cos(lat1_rad) * std::sin(angular_distance) * std::cos(bearing_rad));

            double lon2_rad =
                lon1_rad + std::atan2(std::sin(bearing_rad) * std::sin(angular_distance) * std::cos(lat1_rad),
                                      std::cos(angular_distance) - std::sin(lat1_rad) * std::sin(lat2_rad));

            return WGS{lat2_rad * 180.0 / M_PI, lon2_rad * 180.0 / M_PI, origin.alt};
        }

        // Coordinate validation
        inline bool isValidLatitude(double lat) { return lat >= -90.0 && lat <= 90.0; }

        inline bool isValidLongitude(double lon) { return lon >= -180.0 && lon <= 180.0; }

        inline bool isValidWGS(const WGS &coord) { return isValidLatitude(coord.lat) && isValidLongitude(coord.lon); }

        // Coordinate formatting
        inline std::string formatDMS(double decimal_degrees) {
            bool negative = decimal_degrees < 0;
            decimal_degrees = std::abs(decimal_degrees);

            int degrees = static_cast<int>(decimal_degrees);
            double remainder = (decimal_degrees - degrees) * 60.0;
            int minutes = static_cast<int>(remainder);
            double seconds = (remainder - minutes) * 60.0;

            return (negative ? "-" : "") + std::to_string(degrees) + "°" + std::to_string(minutes) + "'" +
                   std::to_string(seconds) + "\"";
        }

    } // namespace coords

} // namespace concord

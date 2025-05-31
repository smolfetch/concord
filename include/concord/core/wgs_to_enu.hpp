#pragma once

#include "constants.hpp"
#include <tuple>
#include <cmath>

namespace concord {
    
    // Function to convert GPS (lat, lon, alt) to ECEF coordinates
    inline std::tuple<double, double, double> gps_to_ecef(double latitude, double longitude, double altitude) {
        // Convert latitude and longitude to radians
        double cosLat = std::cos(latitude * M_PI / 180.0);
        double sinLat = std::sin(latitude * M_PI / 180.0);
        double cosLong = std::cos(longitude * M_PI / 180.0);
        double sinLong = std::sin(longitude * M_PI / 180.0);

        // Calculate the radius of curvature in the prime vertical
        double N = constants::WGS84_EQUATORIAL_RADIUS / std::sqrt(1 - constants::WGS84_ECCENTRICITY_SQUARED * sinLat * sinLat);

        // Calculate ECEF coordinates
        double x = (N + altitude) * cosLat * cosLong;
        double y = (N + altitude) * cosLat * sinLong;
        double z = (N * (1 - constants::WGS84_ECCENTRICITY_SQUARED) + altitude) * sinLat;

        return std::make_tuple(x, y, z);
    }

    // Function to convert ECEF coordinates to ENU (East, North, Up) with respect to a datum
    inline std::tuple<double, double, double> ecef_to_enu(std::tuple<double, double, double> ecef, std::tuple<double, double, double> datum) {
        double x, y, z;
        std::tie(x, y, z) = ecef;
        double latRef, longRef, altRef;
        std::tie(latRef, longRef, altRef) = datum;

        // Convert the reference latitude and longitude to radians
        double cosLatRef = std::cos(latRef * M_PI / 180.0);
        double sinLatRef = std::sin(latRef * M_PI / 180.0);
        double cosLongRef = std::cos(longRef * M_PI / 180.0);
        double sinLongRef = std::sin(longRef * M_PI / 180.0);

        // Prime vertical radius of curvature at the reference point
        double NRef = constants::WGS84_EQUATORIAL_RADIUS / std::sqrt(1.0 - constants::WGS84_ECCENTRICITY_SQUARED * sinLatRef * sinLatRef);

        // Calculate the reference ECEF coordinates (datum)
        double x0 = (NRef + altRef) * cosLatRef * cosLongRef;
        double y0 = (NRef + altRef) * cosLatRef * sinLongRef;
        double z0 = (NRef * (1 - constants::WGS84_ECCENTRICITY_SQUARED) + altRef) * sinLatRef;

        // Calculate the differences between the ECEF coordinates and the reference ECEF coordinates
        double dx = x - x0;
        double dy = y - y0;
        double dz = z - z0;

        // Calculate the ENU coordinates
        double xEast = -sinLongRef * dx + cosLongRef * dy;
        double yNorth = -cosLongRef * sinLatRef * dx - sinLatRef * sinLongRef * dy + cosLatRef * dz;
        double zUp = cosLatRef * cosLongRef * dx + cosLatRef * sinLongRef * dy + sinLatRef * dz;

        return std::make_tuple(xEast, yNorth, zUp);
    }

    // Function to convert GPS (lat, lon, alt) to ENU (East, North, Up) with respect to a datum
    inline std::tuple<double, double, double> gps_to_enu(double latitude, double longitude, double altitude, double latRef, double longRef, double altRef) {
        // Convert GPS coordinates to ECEF
        std::tuple<double, double, double> ecef = gps_to_ecef(latitude, longitude, altitude);
        // Return the ENU coordinates
        return ecef_to_enu(ecef, std::make_tuple(latRef, longRef, altRef));
    }

    inline std::tuple<double, double, double> enu_to_ecef(std::tuple<double, double, double> enu, std::tuple<double, double, double> datum) {
        // Extract ENU and datum coordinates
        double xEast, yNorth, zUp;
        std::tie(xEast, yNorth, zUp) = enu;
        double latRef, longRef, altRef;
        std::tie(latRef, longRef, altRef) = datum;

        // Compute trigonometric values for the reference latitude and longitude
        double cosLatRef = std::cos(latRef * M_PI / 180);
        double sinLatRef = std::sin(latRef * M_PI / 180);
        double cosLongRef = std::cos(longRef * M_PI / 180);
        double sinLongRef = std::sin(longRef * M_PI / 180);

        // Compute reference ECEF coordinates for the datum
        double cRef = 1 / std::sqrt(cosLatRef * cosLatRef + (1 - constants::WGS84_FLATTENING) * (1 - constants::WGS84_FLATTENING) * sinLatRef * sinLatRef);
        double x0 = (constants::WGS84_EQUATORIAL_RADIUS * cRef + altRef) * cosLatRef * cosLongRef;
        double y0 = (constants::WGS84_EQUATORIAL_RADIUS * cRef + altRef) * cosLatRef * sinLongRef;
        double z0 = (constants::WGS84_EQUATORIAL_RADIUS * cRef * (1 - constants::WGS84_ECCENTRICITY_SQUARED) + altRef) * sinLatRef;

        // Reverse the ENU to ECEF transformation
        double x = x0 + (-sinLongRef * xEast) - (sinLatRef * cosLongRef * yNorth) + (cosLatRef * cosLongRef * zUp);
        double y = y0 + (cosLongRef * xEast) - (sinLatRef * sinLongRef * yNorth) + (cosLatRef * sinLongRef * zUp);
        double z = z0 + (cosLatRef * yNorth) + (sinLatRef * zUp);

        return std::make_tuple(x, y, z);
    }

    inline std::tuple<double, double, double> ecef_to_gps(double x, double y, double z) {
        const double eps = 1e-12; // Convergence threshold
        double longitude = std::atan2(y, x) * 180 / M_PI;

        // Compute initial latitude and height guesses
        double p = std::sqrt(x * x + y * y);
        double latitude = std::atan2(z, p * (1 - constants::WGS84_ECCENTRICITY_SQUARED));
        double N = constants::WGS84_EQUATORIAL_RADIUS / std::sqrt(1 - constants::WGS84_ECCENTRICITY_SQUARED * std::sin(latitude) * std::sin(latitude)); // Prime vertical radius of curvature
        double altitude = p / std::cos(latitude) - N;
        double latOld;

        // Iterate to refine the latitude and height
        do {
            latOld = latitude;
            N = constants::WGS84_EQUATORIAL_RADIUS / std::sqrt(1 - constants::WGS84_ECCENTRICITY_SQUARED * std::sin(latitude) * std::sin(latitude));
            altitude = p / std::cos(latitude) - N;
            latitude = std::atan2(z, p * (1 - constants::WGS84_ECCENTRICITY_SQUARED * N / (N + altitude)));
        } while (std::abs(latitude - latOld) > eps);

        return std::make_tuple(latitude * 180 / M_PI, longitude, altitude);
    }

    inline std::tuple<double, double, double> enu_to_gps(double xEast, double yNorth, double zUp, double latRef, double longRef, double altRef) {
        // Convert ENU to ECEF
        std::tuple<double, double, double> ecef = enu_to_ecef(std::make_tuple(xEast, yNorth, zUp), std::make_tuple(latRef, longRef, altRef));
        // Convert ECEF to GPS
        return ecef_to_gps(std::get<0>(ecef), std::get<1>(ecef), std::get<2>(ecef));
    }
}

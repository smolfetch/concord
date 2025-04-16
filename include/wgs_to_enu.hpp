#pragma once

#include <tuple>
#include <cmath>

// Constants
const double R = 6378137.0;                // Earth's radius in meters (equatorial radius)
const double f = 1.0 / 298.257223563;      // Flattening factor
const double e2 = 2 * f - f * f;           // Square of eccentricity


namespace concord {
    // Function to convert GPS (lat, lon, alt) to ECEF coordinates
    std::tuple<double, double, double> inline gps_to_ecef(double latitude, double longitude, double altitude) {
        // Convert latitude and longitude to radians
        double cosLat = std::cos(latitude * M_PI / 180.0);
        double sinLat = std::sin(latitude * M_PI / 180.0);
        double cosLong = std::cos(longitude * M_PI / 180.0);
        double sinLong = std::sin(longitude * M_PI / 180.0);
        // Prime vertical radius of curvature
        double N = R / std::sqrt(1.0 - e2 * sinLat * sinLat);
        // Calculate ECEF coordinates
        double x = (N + altitude) * cosLat * cosLong;
        double y = (N + altitude) * cosLat * sinLong;
        double z = (N * (1 - e2) + altitude) * sinLat;
        // Return the ECEF coordinates
        return std::make_tuple(x, y, z);
    }

    // Function to convert ECEF coordinates to ENU (East, North, Up) with respect to a datum
    std::tuple<double, double, double> inline ecef_to_enu(std::tuple<double, double, double> ecef, std::tuple<double, double, double> datum) {
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
        double NRef = R / std::sqrt(1.0 - e2 * sinLatRef * sinLatRef);
        // Calculate the reference ECEF coordinates (datum)
        double x0 = (NRef + altRef) * cosLatRef * cosLongRef;
        double y0 = (NRef + altRef) * cosLatRef * sinLongRef;
        double z0 = (NRef * (1 - e2) + altRef) * sinLatRef;
        // Calculate the differences between the ECEF coordinates and the reference ECEF coordinates
        double dx = x - x0;
        double dy = y - y0;
        double dz = z - z0;
        // Calculate the ENU coordinates
        double xEast = -sinLongRef * dx + cosLongRef * dy;
        double yNorth = -cosLongRef * sinLatRef * dx - sinLatRef * sinLongRef * dy + cosLatRef * dz;
        double zUp = cosLatRef * cosLongRef * dx + cosLatRef * sinLongRef * dy + sinLatRef * dz;
        // Return the ENU coordinates
        return std::make_tuple(xEast, yNorth, zUp);
    }

    // Function to convert GPS (lat, lon, alt) to ENU (East, North, Up) with respect to a datum
    std::tuple<double, double, double> inline gps_to_enu(double latitude, double longitude, double altitude, double latRef, double longRef, double altRef) {
        // Convert GPS coordinates to ECEF
        std::tuple<double, double, double> ecef = gps_to_ecef(latitude, longitude, altitude);
        // Return the ENU coordinates
        return ecef_to_enu(ecef, std::make_tuple(latRef, longRef, altRef));
    }

    std::tuple<double, double, double> inline enu_to_ecef(std::tuple<double, double, double> enu, std::tuple<double, double, double> datum) {
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
        double cRef = 1 / std::sqrt(cosLatRef * cosLatRef + (1 - f) * (1 - f) * sinLatRef * sinLatRef);
        double x0 = (R * cRef + altRef) * cosLatRef * cosLongRef;
        double y0 = (R * cRef + altRef) * cosLatRef * sinLongRef;
        double z0 = (R * cRef * (1 - e2) + altRef) * sinLatRef;
        // Reverse the ENU to ECEF transformation
        double x = x0 + (-sinLongRef * xEast) - (sinLatRef * cosLongRef * yNorth) + (cosLatRef * cosLongRef * zUp);
        double y = y0 + (cosLongRef * xEast) - (sinLatRef * sinLongRef * yNorth) + (cosLatRef * sinLongRef * zUp);
        double z = z0 + (cosLatRef * yNorth) + (sinLatRef * zUp);
        // Return the ECEF coordinates
        return std::make_tuple(x, y, z);
    }

    std::tuple<double, double, double> inline ecef_to_gps(double x, double y, double z) {
        const double e2 = f * (2 - f); // Square of eccentricity
        const double eps = 1e-12; // Convergence threshold
        double longitude = std::atan2(y, x) * 180 / M_PI;
        // Compute initial latitude and height guesses
        double p = std::sqrt(x * x + y * y);
        double latitude = std::atan2(z, p * (1 - e2));
        double N = R / std::sqrt(1 - e2 * std::sin(latitude) * std::sin(latitude)); // Prime vertical radius of curvature
        double altitude = p / std::cos(latitude) - N;
        double latOld;
        // Iterate to refine the latitude and height
        do {
            latOld = latitude;
            N = R / std::sqrt(1 - e2 * std::sin(latitude) * std::sin(latitude));
            altitude = p / std::cos(latitude) - N;
            latitude = std::atan2(z, p * (1 - e2 * N / (N + altitude)));
        } while (std::abs(latitude - latOld) > eps);
        // Return the GPS coordinates
        return std::make_tuple(latitude * 180 / M_PI, longitude, altitude);
    }

    std::tuple<double, double, double> inline enu_to_gps(double xEast, double yNorth, double zUp, double latRef, double longRef, double altRef) {
        // Convert ENU to ECEF
        std::tuple<double, double, double> ecef = enu_to_ecef(std::make_tuple(xEast, yNorth, zUp), std::make_tuple(latRef, longRef, altRef));
        // Convert ECEF to GPS
        return ecef_to_gps(std::get<0>(ecef), std::get<1>(ecef), std::get<2>(ecef));
    }
}

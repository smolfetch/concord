#pragma once

#include "types_basic.hpp"
#include "wgs_to_utm.hpp"
#include "wgs_to_enu.hpp"
#include <tuple>
#include <string>
#include <optional>

namespace concord {

    // Extended coordinate system support
    
    // UTM coordinate type
    struct UTM {
        double easting = 0.0;
        double northing = 0.0;
        double altitude = 0.0;
        int zone = 0;
        bool is_north = true;
        
        UTM() = default;
        UTM(double e, double n, double alt, int z, bool north)
            : easting(e), northing(n), altitude(alt), zone(z), is_north(north) {}
        
        // WGS toWGS() const {
        //     auto [lat, lon] = utm_to_wgs(easting, northing, zone, is_north);
        //     return WGS{lat, lon, altitude};
        // }
        
        // static UTM fromWGS(const WGS& wgs) {
        //     auto [e, n, z, north] = wgs_to_utm(wgs.lat, wgs.lon);
        //     return UTM{e, n, wgs.alt, z, north};
        // }
        
        bool is_set() const { return easting != 0.0 || northing != 0.0; }
        
        std::string toString() const {
            return std::to_string(zone) + (is_north ? "N " : "S ") +
                   std::to_string(easting) + " " + std::to_string(northing);
        }
    };
    
    // ECEF (Earth-Centered, Earth-Fixed) coordinate type
    struct ECEF {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        
        ECEF() = default;
        ECEF(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
        
        WGS toWGS() const {
            auto [lat, lon, alt] = ecef_to_gps(x, y, z);
            return WGS{lat, lon, alt};
        }
        
        static ECEF fromWGS(const WGS& wgs) {
            auto [x, y, z] = gps_to_ecef(wgs.lat, wgs.lon, wgs.alt);
            return ECEF{x, y, z};
        }
        
        bool is_set() const { return x != 0.0 || y != 0.0 || z != 0.0; }
    };
    
    // Local Tangent Plane (LTP) coordinate system
    struct LTP {
        double north = 0.0;  // Same as ENU.y
        double east = 0.0;   // Same as ENU.x
        double up = 0.0;     // Same as ENU.z
        
        LTP() = default;
        LTP(double n, double e, double u) : north(n), east(e), up(u) {}
        
        ENU toENU() const {
            return ENU{east, north, up};
        }
        
        static LTP fromENU(const ENU& enu) {
            return LTP{enu.y, enu.x, enu.z};
        }
        
        bool is_set() const { return north != 0.0 || east != 0.0; }
    };
    
    // State Plane Coordinate System (US-specific)
    struct StatePlane {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        int zone = 0;       // State plane zone
        std::string units = "meters";  // "meters" or "feet"
        
        StatePlane() = default;
        StatePlane(double x_, double y_, double z_, int zone_, const std::string& u = "meters")
            : x(x_), y(y_), z(z_), zone(zone_), units(u) {}
        
        bool is_set() const { return x != 0.0 || y != 0.0; }
        
        // Note: Full State Plane conversion would require additional libraries
        // This is a placeholder for the structure
    };
    
    // British National Grid (UK-specific)
    struct BNG {
        double easting = 0.0;
        double northing = 0.0;
        double altitude = 0.0;
        
        BNG() = default;
        BNG(double e, double n, double alt = 0.0) : easting(e), northing(n), altitude(alt) {}
        
        bool is_set() const { return easting != 0.0 || northing != 0.0; }
        
        std::string toGridRef() const {
            // Simplified OS Grid Reference conversion
            // Full implementation would need proper OSGB36 conversion
            int e_100km = static_cast<int>(easting / 100000);
            int n_100km = static_cast<int>(northing / 100000);
            
            // Simplified grid square lookup
            char first_letter = 'S' + (n_100km * 5 + e_100km) / 25;
            char second_letter = 'A' + (n_100km * 5 + e_100km) % 25;
            
            int e_remainder = static_cast<int>(easting) % 100000;
            int n_remainder = static_cast<int>(northing) % 100000;
            
            return std::string(1, first_letter) + std::string(1, second_letter) +
                   std::to_string(e_remainder) + std::to_string(n_remainder);
        }
    };
    
    // Geographic coordinate with different datums
    enum class DatumType {
        WGS84,      // World Geodetic System 1984
        NAD83,      // North American Datum 1983
        NAD27,      // North American Datum 1927
        OSGB36,     // Ordnance Survey Great Britain 1936
        ED50,       // European Datum 1950
        GDA94,      // Geocentric Datum of Australia 1994
        TOKYO       // Tokyo Datum
    };
    
    struct GeographicCoord {
        double latitude = 0.0;
        double longitude = 0.0;
        double altitude = 0.0;
        DatumType datum = DatumType::WGS84;
        
        GeographicCoord() = default;
        GeographicCoord(double lat, double lon, double alt = 0.0, DatumType d = DatumType::WGS84)
            : latitude(lat), longitude(lon), altitude(alt), datum(d) {}
        
        WGS toWGS() const {
            // For now, assume WGS84 or simple conversion
            // Full implementation would need datum transformation parameters
            return WGS{latitude, longitude, altitude};
        }
        
        bool is_set() const { return latitude != 0.0 || longitude != 0.0; }
        
        std::string datumString() const {
            switch (datum) {
                case DatumType::WGS84: return "WGS84";
                case DatumType::NAD83: return "NAD83";
                case DatumType::NAD27: return "NAD27";
                case DatumType::OSGB36: return "OSGB36";
                case DatumType::ED50: return "ED50";
                case DatumType::GDA94: return "GDA94";
                case DatumType::TOKYO: return "TOKYO";
                default: return "UNKNOWN";
            }
        }
    };
    
    // Coordinate system conversion utilities
    namespace coords {
        
        // Distance calculations accounting for Earth's curvature
        double haversineDistance(const WGS& a, const WGS& b) {
            constexpr double R = 6371000.0; // Earth's radius in meters
            
            double lat1_rad = a.lat * M_PI / 180.0;
            double lat2_rad = b.lat * M_PI / 180.0;
            double dlat_rad = (b.lat - a.lat) * M_PI / 180.0;
            double dlon_rad = (b.lon - a.lon) * M_PI / 180.0;
            
            double a_val = std::sin(dlat_rad/2) * std::sin(dlat_rad/2) +
                          std::cos(lat1_rad) * std::cos(lat2_rad) *
                          std::sin(dlon_rad/2) * std::sin(dlon_rad/2);
            
            double c = 2 * std::atan2(std::sqrt(a_val), std::sqrt(1-a_val));
            
            return R * c;
        }
        
        // Bearing calculation
        double bearing(const WGS& from, const WGS& to) {
            double lat1_rad = from.lat * M_PI / 180.0;
            double lat2_rad = to.lat * M_PI / 180.0;
            double dlon_rad = (to.lon - from.lon) * M_PI / 180.0;
            
            double y = std::sin(dlon_rad) * std::cos(lat2_rad);
            double x = std::cos(lat1_rad) * std::sin(lat2_rad) -
                      std::sin(lat1_rad) * std::cos(lat2_rad) * std::cos(dlon_rad);
            
            double bearing_rad = std::atan2(y, x);
            return std::fmod(bearing_rad * 180.0 / M_PI + 360.0, 360.0);
        }
        
        // Point at distance and bearing
        WGS pointAtDistanceAndBearing(const WGS& origin, double distance, double bearing_deg) {
            constexpr double R = 6371000.0; // Earth's radius in meters
            
            double lat1_rad = origin.lat * M_PI / 180.0;
            double lon1_rad = origin.lon * M_PI / 180.0;
            double bearing_rad = bearing_deg * M_PI / 180.0;
            double angular_distance = distance / R;
            
            double lat2_rad = std::asin(std::sin(lat1_rad) * std::cos(angular_distance) +
                                       std::cos(lat1_rad) * std::sin(angular_distance) * std::cos(bearing_rad));
            
            double lon2_rad = lon1_rad + std::atan2(std::sin(bearing_rad) * std::sin(angular_distance) * std::cos(lat1_rad),
                                                   std::cos(angular_distance) - std::sin(lat1_rad) * std::sin(lat2_rad));
            
            return WGS{lat2_rad * 180.0 / M_PI, lon2_rad * 180.0 / M_PI, origin.alt};
        }
        
        // Coordinate validation
        bool isValidLatitude(double lat) {
            return lat >= -90.0 && lat <= 90.0;
        }
        
        bool isValidLongitude(double lon) {
            return lon >= -180.0 && lon <= 180.0;
        }
        
        bool isValidWGS(const WGS& coord) {
            return isValidLatitude(coord.lat) && isValidLongitude(coord.lon);
        }
        
        // Coordinate formatting
        std::string formatDMS(double decimal_degrees) {
            bool negative = decimal_degrees < 0;
            decimal_degrees = std::abs(decimal_degrees);
            
            int degrees = static_cast<int>(decimal_degrees);
            double remainder = (decimal_degrees - degrees) * 60.0;
            int minutes = static_cast<int>(remainder);
            double seconds = (remainder - minutes) * 60.0;
            
            return (negative ? "-" : "") + std::to_string(degrees) + "°" +
                   std::to_string(minutes) + "'" + std::to_string(seconds) + "\"";
        }
        
    } // namespace coords

} // namespace concord

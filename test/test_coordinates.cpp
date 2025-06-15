#include <concord/geographic/coordinate_utils.hpp>
#include <concord/geographic/wgs_to_enu.hpp>
#include <concord/geographic/wgs_to_utm.hpp>
#include <doctest/doctest.h>

using namespace concord;

TEST_CASE("WGS84 to UTM conversion") {
    SUBCASE("Basic coordinate conversion") {
        // Test conversion for a known location (approximate)
        double lat = 40.7128;  // New York City latitude
        double lon = -74.0060; // New York City longitude

        auto [utm_x, utm_y, utm_zone, is_northern] = wgs_to_utm(lat, lon);

        // Basic sanity checks - UTM coordinates should be reasonable
        CHECK(utm_x > 500000); // Typical UTM easting range
        CHECK(utm_x < 800000);
        CHECK(utm_y > 4000000); // Typical UTM northing for this latitude
        CHECK(utm_y < 5000000);
        CHECK(utm_zone > 0);
        CHECK(utm_zone <= 60);
        CHECK(is_northern == true); // NYC is in northern hemisphere
    }

    SUBCASE("UTM zone calculation") {
        // Test various longitudes to check zone calculation
        CHECK(get_utm_zone(0.0) == 31);   // Prime meridian
        CHECK(get_utm_zone(-180.0) == 1); // Date line west
        CHECK(get_utm_zone(174.0) == 60); // Should be zone 60 (close to but before 180)
        CHECK(get_utm_zone(-74.0) == 18); // New York
        CHECK(get_utm_zone(2.0) == 31);   // Paris
    }
}

TEST_CASE("WGS84 to ENU conversion") {
    SUBCASE("Basic ENU conversion") {
        // Reference point (origin for ENU)
        double ref_lat = 40.0;
        double ref_lon = -74.0;
        double ref_alt = 0.0;

        // Target point (same as reference)
        double lat = 40.0;
        double lon = -74.0;
        double alt = 0.0;

        auto [enu_x, enu_y, enu_z] = gps_to_enu(lat, lon, alt, ref_lat, ref_lon, ref_alt);

        // Should be at origin
        CHECK(enu_x == doctest::Approx(0.0).epsilon(0.01));
        CHECK(enu_y == doctest::Approx(0.0).epsilon(0.01));
        CHECK(enu_z == doctest::Approx(0.0).epsilon(0.01));
    }

    SUBCASE("ENU displacement") {
        // Reference point
        double ref_lat = 40.0;
        double ref_lon = -74.0;
        double ref_alt = 0.0;

        // Point slightly north
        double lat = 40.001; // About 111 meters north
        double lon = -74.0;
        double alt = 0.0;

        auto [enu_x, enu_y, enu_z] = gps_to_enu(lat, lon, alt, ref_lat, ref_lon, ref_alt);

        // Should be positive in North (Y) direction
        CHECK(enu_y > 100.0); // Approximately 111 meters
        CHECK(enu_y < 120.0);
        CHECK(std::abs(enu_x) < 1.0); // Should be minimal east displacement
    }
}

TEST_CASE("Coordinate utilities") {
    SUBCASE("ENU distance calculations") {
        ENU p1(0.0, 0.0, 0.0);
        ENU p2(3.0, 4.0, 0.0);

        double dist = p1.distance_to(p2);
        CHECK(dist == doctest::Approx(5.0)); // 3-4-5 triangle
    }

    SUBCASE("Geodetic distance") {
        // Test distance between two well-known points
        WGS nyc(40.7128, -74.0060, 0.0); // NYC
        WGS la(34.0522, -118.2437, 0.0); // LA

        double dist = coords::haversineDistance(nyc, la);

        // Distance NYC to LA is approximately 3944 km
        CHECK(dist > 3900000.0); // meters
        CHECK(dist < 4000000.0);
    }
}

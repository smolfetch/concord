#include <cassert>
#include <cmath>
#include <concord/geographic/wgs_to_enu.hpp>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace concord;

// Test RTK-level precision with high-precision coordinates (8+ decimal places)
void test_rtk_precision() {
    std::cout << "\n=== RTK Precision Test (Sub-Centimeter Level) ===" << std::endl;
    std::cout << std::fixed << std::setprecision(15);

    // Test points with very high precision (8+ decimal places for RTK applications)
    std::vector<std::tuple<double, double, double>> test_points = {
        // Silicon Valley - Google HQ area with sub-cm precision
        {37.42200000000000, -122.08400000000000, 100.000}, // Base reference
        {37.42200001234567, -122.08400001234567, 100.123}, // ~13cm displacement
        {37.42200000012345, -122.08400000012345, 100.012}, // ~1.3cm displacement
        {37.42200000001234, -122.08400000001234, 100.001}, // ~1.3mm displacement

        // London area with high precision
        {51.50740000000000, -0.12780000000000, 50.000}, // Base reference
        {51.50740000123456, -0.12780000123456, 50.012}, // ~1.2cm displacement
        {51.50740000001234, -0.12780000001234, 50.001}, // ~1.2mm displacement

        // Equator/Prime Meridian intersection
        {0.00000000000000, 0.00000000000000, 0.000}, // Origin
        {0.00000001234567, 0.00000001234567, 0.012}, // ~13mm displacement

        // High latitude test (near Arctic)
        {75.00000000000000, 10.00000000000000, 500.000}, // Base reference
        {75.00000001234567, 10.00000001234567, 500.012}  // Small displacement
    };

    double max_position_error_mm = 0.0;
    double max_altitude_error_mm = 0.0;
    int test_count = 0;

    for (size_t i = 0; i < test_points.size(); i += 2) {
        const auto &[lat1, lon1, alt1] = test_points[i];
        const auto &[lat2, lon2, alt2] = test_points[i + 1];

        test_count++;
        std::cout << "\nTest " << test_count << ":" << std::endl;
        std::cout << "  Reference: " << lat1 << ", " << lon1 << ", " << alt1 << std::endl;
        std::cout << "  Target:    " << lat2 << ", " << lon2 << ", " << alt2 << std::endl;

        // Calculate expected displacement
        double lat_diff_m = (lat2 - lat1) * 111000.0; // ~111km per degree latitude
        double lon_diff_m = (lon2 - lon1) * 111000.0 * std::cos(lat1 * M_PI / 180.0);
        double alt_diff_m = alt2 - alt1;
        double expected_distance =
            std::sqrt(lat_diff_m * lat_diff_m + lon_diff_m * lon_diff_m + alt_diff_m * alt_diff_m);

        // Forward conversion: GPS -> ENU
        auto [x, y, z] = gps_to_enu(lat2, lon2, alt2, lat1, lon1, alt1);
        double enu_distance = std::sqrt(x * x + y * y + z * z);

        std::cout << "  Expected distance: " << expected_distance << " m" << std::endl;
        std::cout << "  ENU displacement:  " << x << ", " << y << ", " << z << " m" << std::endl;
        std::cout << "  ENU distance:      " << enu_distance << " m" << std::endl;

        // Reverse conversion: ENU -> GPS (round-trip test)
        auto [lat_back, lon_back, alt_back] = enu_to_gps(x, y, z, lat1, lon1, alt1);

        // Calculate errors
        double lat_error = std::abs(lat2 - lat_back);
        double lon_error = std::abs(lon2 - lon_back);
        double alt_error = std::abs(alt2 - alt_back);

        // Convert coordinate errors to millimeters
        double lat_error_mm = lat_error * 111000.0 * 1000.0; // degrees to mm
        double lon_error_mm = lon_error * 111000.0 * std::cos(lat2 * M_PI / 180.0) * 1000.0;
        double alt_error_mm = alt_error * 1000.0; // meters to mm
        double position_error_mm = std::sqrt(lat_error_mm * lat_error_mm + lon_error_mm * lon_error_mm);

        std::cout << "  Round-trip errors:" << std::endl;
        std::cout << "    Latitude:  " << lat_error << " deg (" << lat_error_mm << " mm)" << std::endl;
        std::cout << "    Longitude: " << lon_error << " deg (" << lon_error_mm << " mm)" << std::endl;
        std::cout << "    Altitude:  " << alt_error << " m (" << alt_error_mm << " mm)" << std::endl;
        std::cout << "    Position:  " << position_error_mm << " mm" << std::endl;

        max_position_error_mm = std::max(max_position_error_mm, position_error_mm);
        max_altitude_error_mm = std::max(max_altitude_error_mm, alt_error_mm);
    }

    std::cout << "\n=== Precision Summary ===" << std::endl;
    std::cout << "Maximum position error: " << max_position_error_mm << " mm" << std::endl;
    std::cout << "Maximum altitude error: " << max_altitude_error_mm << " mm" << std::endl;

    // RTK accuracy requirements: < 10mm horizontal, < 20mm vertical
    bool rtk_horizontal = max_position_error_mm < 10.0; // < 1cm
    bool rtk_vertical = max_altitude_error_mm < 20.0;   // < 2cm

    std::cout << "\nRTK Accuracy Assessment:" << std::endl;
    std::cout << "  Horizontal (< 10mm): " << (rtk_horizontal ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Vertical (< 20mm):   " << (rtk_vertical ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Overall RTK grade:   " << (rtk_horizontal && rtk_vertical ? "PASS" : "FAIL") << std::endl;
}

// Test coordinate transformation consistency
void test_transformation_consistency() {
    std::cout << "\n=== Transformation Consistency Test ===" << std::endl;
    std::cout << std::fixed << std::setprecision(15);

    double lat = 37.422000123456789, lon = -122.084000123456789, alt = 100.123456789;
    double lat_ref = 37.421000123456789, lon_ref = -122.083000123456789, alt_ref = 50.123456789;

    std::cout << "Test coordinates: " << lat << ", " << lon << ", " << alt << std::endl;
    std::cout << "Reference datum:  " << lat_ref << ", " << lon_ref << ", " << alt_ref << std::endl;

    // Path 1: GPS -> ENU directly
    auto [x1, y1, z1] = gps_to_enu(lat, lon, alt, lat_ref, lon_ref, alt_ref);

    // Path 2: GPS -> ECEF -> ENU
    auto [x_ecef, y_ecef, z_ecef] = gps_to_ecef(lat, lon, alt);
    auto [x2, y2, z2] =
        ecef_to_enu(std::make_tuple(x_ecef, y_ecef, z_ecef), std::make_tuple(lat_ref, lon_ref, alt_ref));

    double consistency_error = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));

    std::cout << "Direct GPS->ENU:    " << x1 << ", " << y1 << ", " << z1 << std::endl;
    std::cout << "GPS->ECEF->ENU:     " << x2 << ", " << y2 << ", " << z2 << std::endl;
    std::cout << "Consistency error:  " << consistency_error << " m (" << consistency_error * 1000 << " mm)"
              << std::endl;

    if (consistency_error > 1e-12) {
        std::cout << "WARNING: Significant path inconsistency!" << std::endl;
    } else {
        std::cout << "Path consistency:   EXCELLENT" << std::endl;
    }
}

// Test ECEF round-trip accuracy
void test_ecef_roundtrip() {
    std::cout << "\n=== ECEF Round-trip Test ===" << std::endl;
    std::cout << std::fixed << std::setprecision(15);

    std::vector<std::tuple<double, double, double>> test_coords = {
        {37.422000123456789, -122.084000123456789, 100.123456789},
        {51.507400123456789, -0.127800123456789, 50.123456789},
        {0.000000123456789, 0.000000123456789, 0.123456789},
        {-33.868800123456789, 151.209300123456789, 25.123456789}};

    double max_error_mm = 0.0;

    for (const auto &[lat, lon, alt] : test_coords) {
        std::cout << "\nTesting: " << lat << ", " << lon << ", " << alt << std::endl;

        // GPS -> ECEF -> GPS round trip
        auto [x, y, z] = gps_to_ecef(lat, lon, alt);
        auto [lat_back, lon_back, alt_back] = ecef_to_gps(x, y, z);

        double lat_error = std::abs(lat - lat_back);
        double lon_error = std::abs(lon - lon_back);
        double alt_error = std::abs(alt - alt_back);

        double lat_error_mm = lat_error * 111000.0 * 1000.0;
        double lon_error_mm = lon_error * 111000.0 * std::cos(lat * M_PI / 180.0) * 1000.0;
        double alt_error_mm = alt_error * 1000.0;
        double total_error_mm =
            std::sqrt(lat_error_mm * lat_error_mm + lon_error_mm * lon_error_mm + alt_error_mm * alt_error_mm);

        std::cout << "  ECEF: " << x << ", " << y << ", " << z << std::endl;
        std::cout << "  Recovered: " << lat_back << ", " << lon_back << ", " << alt_back << std::endl;
        std::cout << "  Total error: " << total_error_mm << " mm" << std::endl;

        max_error_mm = std::max(max_error_mm, total_error_mm);
    }

    std::cout << "\nMaximum ECEF round-trip error: " << max_error_mm << " mm" << std::endl;
    std::cout << "ECEF precision: " << (max_error_mm < 0.1 ? "EXCELLENT" : "GOOD") << std::endl;
}

// Test small displacement accuracy (critical for RTK applications)
void test_small_displacements() {
    std::cout << "\n=== Small Displacement Test ===" << std::endl;
    std::cout << std::fixed << std::setprecision(15);

    double base_lat = 37.422000000000000;
    double base_lon = -122.084000000000000;
    double base_alt = 100.000000000000000;

    // Test displacements at different scales
    std::vector<std::tuple<std::string, double, double, double>> displacements = {
        {"1mm N", 0.000000009, 0.0, 0.0},  // ~1mm north
        {"1mm E", 0.0, 0.000000013, 0.0},  // ~1mm east
        {"1mm U", 0.0, 0.0, 0.001},        // 1mm up
        {"1cm N", 0.000000090, 0.0, 0.0},  // ~1cm north
        {"1cm E", 0.0, 0.000000130, 0.0},  // ~1cm east
        {"1cm U", 0.0, 0.0, 0.01},         // 1cm up
        {"10cm N", 0.000000900, 0.0, 0.0}, // ~10cm north
        {"10cm E", 0.0, 0.000001300, 0.0}, // ~10cm east
        {"10cm U", 0.0, 0.0, 0.1}          // 10cm up
    };

    for (const auto &[name, dlat, dlon, dalt] : displacements) {
        double target_lat = base_lat + dlat;
        double target_lon = base_lon + dlon;
        double target_alt = base_alt + dalt;

        std::cout << "\nTesting " << name << " displacement:" << std::endl;

        // Convert to ENU
        auto [x, y, z] = gps_to_enu(target_lat, target_lon, target_alt, base_lat, base_lon, base_alt);
        double enu_magnitude = std::sqrt(x * x + y * y + z * z);

        // Convert back to GPS
        auto [lat_back, lon_back, alt_back] = enu_to_gps(x, y, z, base_lat, base_lon, base_alt);

        // Calculate errors
        double lat_error_mm = std::abs(target_lat - lat_back) * 111000.0 * 1000.0;
        double lon_error_mm = std::abs(target_lon - lon_back) * 111000.0 * std::cos(target_lat * M_PI / 180.0) * 1000.0;
        double alt_error_mm = std::abs(target_alt - alt_back) * 1000.0;

        std::cout << "  ENU: " << x * 1000 << ", " << y * 1000 << ", " << z * 1000 << " mm" << std::endl;
        std::cout << "  ENU magnitude: " << enu_magnitude * 1000 << " mm" << std::endl;
        std::cout << "  Round-trip error: " << lat_error_mm << ", " << lon_error_mm << ", " << alt_error_mm << " mm"
                  << std::endl;
    }
}

int main() {
    std::cout << "=== Concord WGS84 ↔ ENU Precision Test Suite ===" << std::endl;
    std::cout << "Testing coordinate conversion precision for RTK-level applications" << std::endl;

    test_rtk_precision();
    test_transformation_consistency();
    test_ecef_roundtrip();
    test_small_displacements();

    std::cout << "\n=== Test Suite Complete ===" << std::endl;
    std::cout << "For RTK applications, ensure all tests show PASS status" << std::endl;
    std::cout << "and errors are well below 10mm horizontal, 20mm vertical" << std::endl;

    return 0;
}

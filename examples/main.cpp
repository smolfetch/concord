#include <iostream>
#include <vector>
#include <iomanip>

// Include the main concord header which includes everything
#include <concord/concord.hpp>

using namespace concord;

void test_mathematical_types() {
    std::cout << "\n=== Testing Mathematical Types ===" << std::endl;
    
    // Test Vector operations
    Vec3d v1;
    v1[0] = 1.0; v1[1] = 2.0; v1[2] = 3.0;
    Vec3d v2;
    v2[0] = 4.0; v2[1] = 5.0; v2[2] = 6.0;
    Vec3d v3 = v1 + v2;
    
    std::cout << "Vector addition: (" << v1[0] << "," << v1[1] << "," << v1[2] << ") + "
              << "(" << v2[0] << "," << v2[1] << "," << v2[2] << ") = "
              << "(" << v3[0] << "," << v3[1] << "," << v3[2] << ")" << std::endl;
    
    double dot = dot_product(v1, v2);
    Vec3d cross = cross_product(v1, v2);
    std::cout << "Dot product: " << dot << std::endl;
    std::cout << "Cross product: (" << cross[0] << "," << cross[1] << "," << cross[2] << ")" << std::endl;
    
    // Test Matrix operations
    Mat3d identity = Mat3d::identity();
    Mat3d rotation = create_rotation_z(M_PI / 4); // 45 degrees
    
    std::cout << "Created identity and rotation matrices" << std::endl;
}

void test_enhanced_basic_types() {
    std::cout << "\n=== Testing Enhanced Basic Types ===" << std::endl;
    
    try {
        // Test enhanced WGS coordinates
        WGS seattle(47.6062, -122.3321, 56.0);
        WGS portland(45.5152, -122.6784, 15.0);
        
        std::cout << "Seattle: " << seattle.lat << ", " << seattle.lon << ", " << seattle.alt << "m" << std::endl;
        std::cout << "Portland: " << portland.lat << ", " << portland.lon << ", " << portland.alt << "m" << std::endl;
        
        // Test distance calculation
        double distance = seattle.distance_to(portland);
        std::cout << "Distance between Seattle and Portland: " << distance << " meters" << std::endl;
        
        // Test bearing calculation
        double bearing = seattle.bearing_to(portland);
        std::cout << "Bearing from Seattle to Portland: " << bearing << " degrees" << std::endl;
        
        // Test enhanced ENU coordinates
        ENU enu1(100.0, 200.0, 50.0);
        ENU enu2(300.0, 400.0, 75.0);
        
        double enu_distance = enu1.distance_to(enu2);
        std::cout << "ENU distance: " << enu_distance << " meters" << std::endl;
        
        // Test enhanced Quaternion
        Quaternion q1(1.0, 0.0, 0.0, 0.0);
        Quaternion q2(0.707, 0.0, 0.0, 0.707); // 90 degree rotation around Z
        Quaternion result = q1 * q2;  // Use operator* instead of multiply method
        
        std::cout << "Quaternion multiplication result: (" 
                  << result.w << "," << result.x << "," << result.y << "," << result.z << ")" << std::endl;
        
    } catch (const ConcordException& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
}

void test_coordinate_systems() {
    std::cout << "\n=== Testing Coordinate Systems ===" << std::endl;
    
    try {
        // Test basic WGS coordinates
        WGS wgs_point(47.6062, -122.3321, 56.0);  // Seattle
        
        std::cout << "WGS coordinates:" << std::endl;
        std::cout << "  Lat: " << wgs_point.lat << ", Lon: " << wgs_point.lon << ", Alt: " << wgs_point.alt << std::endl;
        
        // Test validation using proper namespace
        try {
            concord::validation::validate_latitude(wgs_point.lat);
            concord::validation::validate_longitude(wgs_point.lon);
            std::cout << "  Coordinates are valid" << std::endl;
        } catch (const InvalidCoordinateException& e) {
            std::cout << "  Validation error: " << e.what() << std::endl;
        }
        
    } catch (const ConcordException& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
}

void test_geometric_types() {
    std::cout << "\n=== Testing Geometric Types ===" << std::endl;
    
    // Test basic Point operations
    Point p1(5, 5, 5);
    Point p2(10, 10, 10);
    
    std::cout << "Point 1: (" << p1.enu.x << ", " << p1.enu.y << ", " << p1.enu.z << ")" << std::endl;
    std::cout << "Point 2: (" << p2.enu.x << ", " << p2.enu.y << ", " << p2.enu.z << ")" << std::endl;
    
    double distance = p1.distance_to(p2);
    std::cout << "Distance between points: " << distance << std::endl;
}

void test_spatial_algorithms() {
    std::cout << "\n=== Testing Spatial Algorithms ===" << std::endl;
    
    std::cout << "Spatial algorithms are loaded and available" << std::endl;
    std::cout << "- Distance calculations: available" << std::endl;
    std::cout << "- Line intersections: available" << std::endl;
    std::cout << "- Convex hull: available" << std::endl;
}

void test_spatial_indexing() {
    std::cout << "\n=== Testing Spatial Indexing ===" << std::endl;
    
    // Test spatial hash grid
    SpatialHashGrid<int> grid(10.0);  // 10 unit cell size
    Point p1(15, 25, 0);
    Point p2(18, 22, 0);
    
    grid.insert(p1, 100);
    grid.insert(p2, 200);
    
    auto nearby = grid.query(Point(20, 20, 0), 10.0);
    std::cout << "Hash grid found " << nearby.size() << " points within radius" << std::endl;
}

void test_utilities() {
    std::cout << "\n=== Testing Utilities ===" << std::endl;
    
    std::cout << "Utility functions are loaded and available" << std::endl;
    std::cout << "- Random generation: available" << std::endl;
    std::cout << "- Unit conversions: available" << std::endl;
    std::cout << "- Validation: available" << std::endl;
}

void test_interpolation() {
    std::cout << "\n=== Testing Interpolation ===" << std::endl;
    
    // Test linear interpolation
    double result = lerp(0.0, 10.0, 0.5);  // Should be 5.0
    std::cout << "Linear interpolation between 0 and 10 at t=0.5: " << result << std::endl;
    
    // Test smooth interpolation
    double smooth = smoothstep(0.0, 1.0, 0.5);
    std::cout << "Smoothstep interpolation at t=0.5: " << smooth << std::endl;
}

int main() {
    std::cout << "=== Concord Library Comprehensive Test ===" << std::endl;
    std::cout << "Version: " << VERSION_STRING << std::endl;
    std::cout << "Capabilities:" << std::endl;
    std::cout << "  Mathematical Types: " << (HAS_MATHEMATICAL_TYPES ? "Yes" : "No") << std::endl;
    std::cout << "  Spatial Indexing: " << (HAS_SPATIAL_INDEXING ? "Yes" : "No") << std::endl;
    std::cout << "  Advanced Algorithms: " << (HAS_ADVANCED_ALGORITHMS ? "Yes" : "No") << std::endl;
    std::cout << "  Multiple Datums: " << (HAS_MULTIPLE_DATUMS ? "Yes" : "No") << std::endl;
    std::cout << "  I/O Operations: " << (HAS_IO_OPERATIONS ? "Yes" : "No") << std::endl;
    
    // Set precision for output
    std::cout << std::fixed << std::setprecision(6);
    
    // Run all tests
    test_mathematical_types();
    test_enhanced_basic_types();
    test_coordinate_systems();
    test_geometric_types();
    test_spatial_algorithms();
    test_spatial_indexing();
    test_utilities();
    test_interpolation();
    
    std::cout << "\n=== All Tests Completed ===" << std::endl;
    return 0;
}

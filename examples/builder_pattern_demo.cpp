#include <concord/concord.hpp>
#include <iomanip>
#include <iostream>

using namespace concord;

int main() {
    std::cout << std::fixed << std::setprecision(6);

    // Example 1: Basic Point to ENU to WGS transformation
    std::cout << "=== Builder Pattern Examples ===" << std::endl;

    // Define a reference datum (Seattle area)
    Datum seattle_datum(47.6062, -122.3321, 56.0);

    // Example 1: Point -> ENU -> WGS
    std::cout << "\n1. Point -> ENU -> WGS transformation:" << std::endl;
    Point local_point(100.0, 200.0, 50.0);
    std::cout << "   Local Point: (" << local_point.x << ", " << local_point.y << ", " << local_point.z << ")"
              << std::endl;

    // Using template-based builder pattern with fluent chaining
    auto wgs_result = convert(local_point)
        .withDatum(seattle_datum)
        .as<ENU>()
        .to<WGS>()
        .build();

    std::cout << "   WGS Result: (" << wgs_result.lat << ", " << wgs_result.lon << ", " << wgs_result.alt << ")"
              << std::endl;

    // Example 2: WGS -> ENU -> Point chain
    std::cout << "\n2. WGS -> ENU -> Point transformation:" << std::endl;
    WGS portland(45.5152, -122.6784, 15.0);
    std::cout << "   WGS Point: (" << portland.lat << ", " << portland.lon << ", " << portland.alt << ")" << std::endl;

    auto enu_result = convert(portland)
        .withDatum(seattle_datum)
        .as<ENU>()
        .get();  // Extract the ENU result directly

    std::cout << "   ENU Result: (" << enu_result.x << ", " << enu_result.y << ", " << enu_result.z << ")" << std::endl;    // Example 3: Fluent chaining with multiple transformations
    std::cout << "\n3. Multiple transformations in one chain:" << std::endl;

    // Beautiful template-based fluent syntax:
    auto multi_result = convert(local_point)
        .withDatum(seattle_datum)
        .as<ENU>()
        .to<WGS>()
        .get();
    
    std::cout << "   Chained Point->ENU->WGS: (" << multi_result.lat << ", " 
              << multi_result.lon << ", " << multi_result.alt << ")" << std::endl;

    // Example 4: Distance calculations between different coordinate systems
    std::cout << "\n4. Distance calculations:" << std::endl;
    WGS seattle(47.6062, -122.3321, 56.0);
    WGS vancouver(49.2827, -123.1207, 70.0);    // Convert both to ENU for distance calculation
    auto seattle_enu = convert(seattle)
        .withDatum(seattle_datum)
        .as<ENU>()
        .get();
    
    auto vancouver_enu = convert(vancouver)
        .withDatum(seattle_datum)
        .as<ENU>()
        .get();
    
    double distance = seattle_enu.distance_to(vancouver_enu);
    std::cout << "   Distance Seattle->Vancouver: " << distance << " meters" << std::endl;

    return 0;
}

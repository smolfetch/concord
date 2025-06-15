#include <concord/concord.hpp>
#include <iostream>
#include <iomanip>

using namespace concord;

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "=== Template-Based Builder Pattern Examples ===" << std::endl;
    
    // Define a reference datum (Seattle area)
    Datum seattle_datum(47.6062, -122.3321, 56.0);
    
    // Example 1: Point -> ENU using template syntax
    std::cout << "\n1. Point -> ENU using as<ENU>():" << std::endl;
    Point local_point(100.0, 200.0, 50.0);
    std::cout << "   Local Point: (" << local_point.x << ", " << local_point.y << ", " << local_point.z << ")" << std::endl;
    
    auto enu_result = convert(local_point)
        .withDatum(seattle_datum)
        .as<ENU>()
        .get();
    
    std::cout << "   ENU Result: (" << enu_result.x << ", " << enu_result.y << ", " << enu_result.z << ")" << std::endl;
    
    // Example 2: Point -> ENU -> WGS chain with templates
    std::cout << "\n2. Point -> ENU -> WGS using to<>() chain:" << std::endl;
    
    auto wgs_result = convert(local_point)
        .withDatum(seattle_datum)
        .as<ENU>()
        .to<WGS>()
        .build();
    
    std::cout << "   WGS Result: (" << wgs_result.lat << ", " << wgs_result.lon << ", " << wgs_result.alt << ")" << std::endl;
    
    // Example 3: WGS -> ENU using template syntax
    std::cout << "\n3. WGS -> ENU using to<ENU>():" << std::endl;
    WGS portland(45.5152, -122.6784, 15.0);
    std::cout << "   WGS Point: (" << portland.lat << ", " << portland.lon << ", " << portland.alt << ")" << std::endl;
    
    auto enu_from_wgs = convert(portland)
        .withDatum(seattle_datum)
        .to<ENU>()
        .build();
    
    std::cout << "   ENU Result: (" << enu_from_wgs.x << ", " << enu_from_wgs.y << ", " << enu_from_wgs.z << ")" << std::endl;
    
    // Example 4: Complex transformation chain
    std::cout << "\n4. Complex chain: Point -> ENU -> WGS -> ENU -> Point:" << std::endl;
    
    auto final_point = convert(local_point)
        .withDatum(seattle_datum)
        .as<ENU>()
        .to<WGS>()
        .to<ENU>()
        .as<Point>()
        .get();
    
    std::cout << "   Original Point: (" << local_point.x << ", " << local_point.y << ", " << local_point.z << ")" << std::endl;
    std::cout << "   Final Point:    (" << final_point.x << ", " << final_point.y << ", " << final_point.z << ")" << std::endl;
    
    // Example 5: Different factory functions
    std::cout << "\n5. Alternative factory functions:" << std::endl;
    
    auto coord_result = coordinate(local_point).withDatum(seattle_datum).as<ENU>().get();
    auto transform_result = transform(local_point).withDatum(seattle_datum).as<ENU>().get();
    
    std::cout << "   coordinate() factory: (" << coord_result.x << ", " << coord_result.y << ", " << coord_result.z << ")" << std::endl;
    std::cout << "   transform() factory:  (" << transform_result.x << ", " << transform_result.y << ", " << transform_result.z << ")" << std::endl;
    
    std::cout << "\n✅ All template-based transformations completed successfully!" << std::endl;
    
    return 0;
}

#include <concord/concord.hpp>
#include <iostream>

using namespace concord;

int main() {
    std::cout << "Testing new Concord structure..." << std::endl;

    // Test basic Point creation
    Point p(10.0, 20.0, 30.0);
    std::cout << "Point: (" << p.x << ", " << p.y << ", " << p.z << ")" << std::endl;

    // Test WGS coordinates
    try {
        WGS seattle(47.6062, -122.3321, 56.0);
        std::cout << "WGS Seattle: (" << seattle.lat << ", " << seattle.lon << ", " << seattle.alt << ")" << std::endl;

        // Test datum
        Datum datum(47.6062, -122.3321, 56.0);

        // Test ENU creation
        ENU enu_point(100.0, 200.0, 50.0, datum);
        std::cout << "ENU Point: (" << enu_point.x << ", " << enu_point.y << ", " << enu_point.z << ")" << std::endl;

        std::cout << "All basic tests passed!" << std::endl;

    } catch (const std::exception &e) {
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

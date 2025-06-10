#include <doctest/doctest.h>
#include <concord/geometry/circle.hpp>
#include <concord/geometry/rectangle.hpp>
#include <concord/geometry/line.hpp>

using namespace concord;

TEST_CASE("Circle geometry") {
    SUBCASE("Circle area calculation") {
        Point center;
        center.enu.x = 0.0; center.enu.y = 0.0; center.enu.z = 0.0;
        
        Circle circle(center, 5.0);
        CHECK(circle.area() == doctest::Approx(78.5398).epsilon(0.001)); // π * 5^2
    }
    
    SUBCASE("Circle circumference") {
        Point center;
        center.enu.x = 2.0; center.enu.y = 3.0; center.enu.z = 0.0;
        
        Circle circle(center, 3.0);
        CHECK(circle.circumference() == doctest::Approx(18.8496).epsilon(0.001)); // 2 * π * 3
    }
    
    SUBCASE("Point inside circle") {
        Point center;
        center.enu.x = 0.0; center.enu.y = 0.0; center.enu.z = 0.0;
        
        Point test_point;
        test_point.enu.x = 2.0; test_point.enu.y = 2.0; test_point.enu.z = 0.0;
        
        Circle circle(center, 5.0);
        CHECK(circle.contains(test_point) == true);
        
        test_point.enu.x = 6.0; test_point.enu.y = 0.0; test_point.enu.z = 0.0;
        CHECK(circle.contains(test_point) == false);
    }
}

TEST_CASE("Rectangle geometry") {
    SUBCASE("Rectangle area") {
        Point tl, tr, bl, br;
        tl.enu.x = 1.0; tl.enu.y = 4.0; tl.enu.z = 0.0; // top-left
        tr.enu.x = 5.0; tr.enu.y = 4.0; tr.enu.z = 0.0; // top-right
        bl.enu.x = 1.0; bl.enu.y = 1.0; bl.enu.z = 0.0; // bottom-left
        br.enu.x = 5.0; br.enu.y = 1.0; br.enu.z = 0.0; // bottom-right
        
        Rectangle rect(tl, tr, bl, br);
        CHECK(rect.area() == doctest::Approx(12.0)); // 4 * 3
    }
    
    SUBCASE("Rectangle perimeter") {
        Point tl, tr, bl, br;
        tl.enu.x = 1.0; tl.enu.y = 5.0; tl.enu.z = 0.0; // top-left
        tr.enu.x = 5.0; tr.enu.y = 5.0; tr.enu.z = 0.0; // top-right
        bl.enu.x = 1.0; bl.enu.y = 1.0; bl.enu.z = 0.0; // bottom-left
        br.enu.x = 5.0; br.enu.y = 1.0; br.enu.z = 0.0; // bottom-right
        
        Rectangle rect(tl, tr, bl, br);
        CHECK(rect.perimeter() == doctest::Approx(16.0)); // 2 * (4 + 4)
    }
}

TEST_CASE("Line geometry") {
    SUBCASE("Line length calculation") {
        Point start, end;
        start.enu.x = 0.0; start.enu.y = 0.0; start.enu.z = 0.0;
        end.enu.x = 3.0; end.enu.y = 4.0; end.enu.z = 0.0;
        
        Line line(start, end);
        CHECK(line.length() == doctest::Approx(5.0)); // 3-4-5 triangle
    }
}

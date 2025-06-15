#include <concord/geometry/primitives/circle.hpp>
#include <concord/geometry/primitives/line.hpp>
#include <concord/geometry/primitives/rectangle.hpp>
#include <doctest/doctest.h>

using namespace concord;

TEST_CASE("Circle geometry") {
    SUBCASE("Circle area calculation") {
        Point center;
        center.x = 0.0;
        center.y = 0.0;
        center.z = 0.0;

        Circle circle(center, 5.0);
        CHECK(circle.area() == doctest::Approx(78.5398).epsilon(0.001)); // π * 5^2
    }

    SUBCASE("Circle circumference") {
        Point center;
        center.x = 2.0;
        center.y = 3.0;
        center.z = 0.0;

        Circle circle(center, 3.0);
        CHECK(circle.circumference() == doctest::Approx(18.8496).epsilon(0.001)); // 2 * π * 3
    }

    SUBCASE("Point inside circle") {
        Point center;
        center.x = 0.0;
        center.y = 0.0;
        center.z = 0.0;

        Point test_point;
        test_point.x = 2.0;
        test_point.y = 2.0;
        test_point.z = 0.0;

        Circle circle(center, 5.0);
        CHECK(circle.contains(test_point) == true);

        test_point.x = 6.0;
        test_point.y = 0.0;
        test_point.z = 0.0;
        CHECK(circle.contains(test_point) == false);
    }
}

TEST_CASE("Rectangle geometry") {
    SUBCASE("Rectangle area") {
        Point tl, tr, bl, br;
        tl.x = 1.0;
        tl.y = 4.0;
        tl.z = 0.0; // top-left
        tr.x = 5.0;
        tr.y = 4.0;
        tr.z = 0.0; // top-right
        bl.x = 1.0;
        bl.y = 1.0;
        bl.z = 0.0; // bottom-left
        br.x = 5.0;
        br.y = 1.0;
        br.z = 0.0; // bottom-right

        Rectangle rect(tl, tr, bl, br);
        CHECK(rect.area() == doctest::Approx(12.0)); // 4 * 3
    }

    SUBCASE("Rectangle perimeter") {
        Point tl, tr, bl, br;
        tl.x = 1.0;
        tl.y = 5.0;
        tl.z = 0.0; // top-left
        tr.x = 5.0;
        tr.y = 5.0;
        tr.z = 0.0; // top-right
        bl.x = 1.0;
        bl.y = 1.0;
        bl.z = 0.0; // bottom-left
        br.x = 5.0;
        br.y = 1.0;
        br.z = 0.0; // bottom-right

        Rectangle rect(tl, tr, bl, br);
        CHECK(rect.perimeter() == doctest::Approx(16.0)); // 2 * (4 + 4)
    }
}

TEST_CASE("Line geometry") {
    SUBCASE("Line length calculation") {
        Point start, end;
        start.x = 0.0;
        start.y = 0.0;
        start.z = 0.0;
        end.x = 3.0;
        end.y = 4.0;
        end.z = 0.0;

        Line line(start, end);
        CHECK(line.length() == doctest::Approx(5.0)); // 3-4-5 triangle
    }
}

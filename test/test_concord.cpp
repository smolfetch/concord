#include <concord/concord.hpp>
#include <doctest/doctest.h>

using namespace concord;

TEST_CASE("Basic Vector operations") {
    SUBCASE("Vector creation and access") {
        Vector<double, 3> v{1.0, 2.0, 3.0};
        CHECK(v[0] == 1.0);
        CHECK(v[1] == 2.0);
        CHECK(v[2] == 3.0);
    }

    SUBCASE("Vector addition") {
        Vector<double, 3> v1{1.0, 2.0, 3.0};
        Vector<double, 3> v2{4.0, 5.0, 6.0};
        auto result = v1 + v2;
        CHECK(result[0] == 5.0);
        CHECK(result[1] == 7.0);
        CHECK(result[2] == 9.0);
    }

    SUBCASE("Vector subtraction") {
        Vector<double, 3> v1{4.0, 5.0, 6.0};
        Vector<double, 3> v2{1.0, 2.0, 3.0};
        auto result = v1 - v2;
        CHECK(result[0] == 3.0);
        CHECK(result[1] == 3.0);
        CHECK(result[2] == 3.0);
    }
}

TEST_CASE("Circle geometry tests") {
    SUBCASE("Circle creation and basic properties") {
        Point center;
        center.x = 0.0;
        center.y = 0.0;
        center.z = 0.0;

        Circle circle(center, 5.0);

        CHECK(circle.getCenter().x == 0.0);
        CHECK(circle.getCenter().y == 0.0);
        CHECK(doctest::Approx(circle.area()) == M_PI * 25.0);
        CHECK(doctest::Approx(circle.circumference()) == 2 * M_PI * 5.0);
    }

    SUBCASE("Point containment") {
        Point center;
        center.x = 0.0;
        center.y = 0.0;
        center.z = 0.0;

        Circle circle(center, 5.0);

        Point inside;
        inside.x = 2.0;
        inside.y = 2.0;
        inside.z = 0.0;

        Point outside;
        outside.x = 10.0;
        outside.y = 10.0;
        outside.z = 0.0;

        CHECK(circle.contains(inside) == true);
        CHECK(circle.contains(outside) == false);
    }
}

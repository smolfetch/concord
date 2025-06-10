#include <concord/math/math.hpp>
#include <doctest/doctest.h>

using namespace concord;

TEST_CASE("Vector math operations") {
    SUBCASE("Dot product") {
        Vector<double, 3> v1{1.0, 2.0, 3.0};
        Vector<double, 3> v2{4.0, 5.0, 6.0};
        double result = v1.dot(v2);
        CHECK(result == doctest::Approx(32.0)); // 1*4 + 2*5 + 3*6 = 32
    }

    SUBCASE("Vector magnitude") {
        Vector<double, 3> v{3.0, 4.0, 0.0};
        double mag = v.magnitude();
        CHECK(mag == doctest::Approx(5.0)); // sqrt(9 + 16) = 5
    }

    SUBCASE("Vector normalization") {
        Vector<double, 3> v{3.0, 4.0, 0.0};
        auto normalized = v.normalized();
        CHECK(normalized.magnitude() == doctest::Approx(1.0));
        CHECK(normalized[0] == doctest::Approx(0.6));
        CHECK(normalized[1] == doctest::Approx(0.8));
    }

    SUBCASE("Cross product") {
        Vector<double, 3> v1{1.0, 0.0, 0.0};
        Vector<double, 3> v2{0.0, 1.0, 0.0};
        auto result = cross(v1, v2);
        CHECK(result[0] == doctest::Approx(0.0));
        CHECK(result[1] == doctest::Approx(0.0));
        CHECK(result[2] == doctest::Approx(1.0));
    }
}

TEST_CASE("Matrix operations") {
    SUBCASE("Matrix creation and access") {
        Matrix<double, 2, 2> m;
        m[0][0] = 1.0;
        m[0][1] = 2.0;
        m[1][0] = 3.0;
        m[1][1] = 4.0;
        CHECK(m[0][0] == 1.0);
        CHECK(m[0][1] == 2.0);
        CHECK(m[1][0] == 3.0);
        CHECK(m[1][1] == 4.0);
    }

    SUBCASE("Matrix multiplication") {
        Matrix<double, 2, 2> m1;
        m1[0][0] = 1.0;
        m1[0][1] = 2.0;
        m1[1][0] = 3.0;
        m1[1][1] = 4.0;

        Matrix<double, 2, 2> m2;
        m2[0][0] = 5.0;
        m2[0][1] = 6.0;
        m2[1][0] = 7.0;
        m2[1][1] = 8.0;

        auto result = m1 * m2;
        CHECK(result[0][0] == doctest::Approx(19.0)); // 1*5 + 2*7
        CHECK(result[0][1] == doctest::Approx(22.0)); // 1*6 + 2*8
        CHECK(result[1][0] == doctest::Approx(43.0)); // 3*5 + 4*7
        CHECK(result[1][1] == doctest::Approx(50.0)); // 3*6 + 4*8
    }
}

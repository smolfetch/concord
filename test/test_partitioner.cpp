#include <concord/geometry/polygon/partitioner.hpp>
#include <doctest/doctest.h>
#include <vector>

using namespace concord;

TEST_CASE("Polygon Partitioner") {
    SUBCASE("Basic partitioning by area") {
        // Create an L-shaped polygon
        std::vector<Point> l_shape_points;

        // Points of an L-shape
        Point p;
        p.x = 0.0;
        p.y = 0.0;
        p.z = 0.0;
        l_shape_points.push_back(p);
        p.x = 20.0;
        p.y = 0.0;
        p.z = 0.0;
        l_shape_points.push_back(p);
        p.x = 20.0;
        p.y = 10.0;
        p.z = 0.0;
        l_shape_points.push_back(p);
        p.x = 10.0;
        p.y = 10.0;
        p.z = 0.0;
        l_shape_points.push_back(p);
        p.x = 10.0;
        p.y = 20.0;
        p.z = 0.0;
        l_shape_points.push_back(p);
        p.x = 0.0;
        p.y = 20.0;
        p.z = 0.0;
        l_shape_points.push_back(p);

        Polygon l_shape(l_shape_points);

        // Area of the L-shape should be 300 square units
        CHECK(doctest::Approx(l_shape.area()) == 300.0);

        // Create a partitioner
        Partitioner partitioner(l_shape);

        // Partition with a very small area threshold to ensure partition happens
        auto partitioned = partitioner.partition(100.0);

        // Should generate multiple smaller polygons
        CHECK(partitioned.size() > 1);

        // Check that each partitioned polygon has area less than the threshold
        for (const auto &part : partitioned) {
            CHECK(part.area() <= 100.0);
        }

        // Total area of parts should equal original area (approximately)
        double total_area = 0.0;
        for (const auto &part : partitioned) {
            total_area += part.area();
        }
        CHECK(doctest::Approx(total_area) == l_shape.area());
    }

    SUBCASE("Partitioning by convexity") {
        // Create a concave polygon (U-shape)
        std::vector<Point> u_shape_points;

        Point p;
        p.x = 0.0;
        p.y = 0.0;
        p.z = 0.0;
        u_shape_points.push_back(p);
        p.x = 20.0;
        p.y = 0.0;
        p.z = 0.0;
        u_shape_points.push_back(p);
        p.x = 20.0;
        p.y = 20.0;
        p.z = 0.0;
        u_shape_points.push_back(p);
        p.x = 15.0;
        p.y = 20.0;
        p.z = 0.0;
        u_shape_points.push_back(p);
        p.x = 15.0;
        p.y = 5.0;
        p.z = 0.0;
        u_shape_points.push_back(p);
        p.x = 5.0;
        p.y = 5.0;
        p.z = 0.0;
        u_shape_points.push_back(p);
        p.x = 5.0;
        p.y = 20.0;
        p.z = 0.0;
        u_shape_points.push_back(p);
        p.x = 0.0;
        p.y = 20.0;
        p.z = 0.0;
        u_shape_points.push_back(p);

        Polygon u_shape(u_shape_points);

        // This U-shape should have low convexity
        Partitioner partitioner(u_shape);

        // Set criteria to focus on convexity
        Partitioner::PartitionCriteria criteria;
        criteria.max_area = 1000.0;   // High area threshold so it doesn't partition by area
        criteria.min_convexity = 0.8; // High convexity requirement to force partition

        auto partitioned = partitioner.partition(1000.0, criteria);

        // Should split into multiple parts to improve convexity
        CHECK(partitioned.size() > 1);

        // Check that each part has better convexity than original
        double original_convexity = 0.0;
        {
            Polygon convex_hull = spatial::convexHull(u_shape.getPoints());
            double poly_area = u_shape.area();
            double hull_area = convex_hull.area();

            if (hull_area > 1e-10) {
                original_convexity = poly_area / hull_area;
            }
        }

        for (const auto &part : partitioned) {
            Polygon part_convex_hull = spatial::convexHull(part.getPoints());
            double part_area = part.area();
            double part_hull_area = part_convex_hull.area();

            if (part_hull_area > 1e-10) {
                double part_convexity = part_area / part_hull_area;
                // Each part should have better convexity than original
                CHECK(part_convexity > original_convexity);
            }
        }
    }

    SUBCASE("Aspect ratio partitioning") {
        // Create a very elongated rectangle
        std::vector<Point> rect_points;

        Point p;
        p.x = 0.0;
        p.y = 0.0;
        p.z = 0.0;
        rect_points.push_back(p);
        p.x = 40.0;
        p.y = 0.0;
        p.z = 0.0;
        rect_points.push_back(p);
        p.x = 40.0;
        p.y = 5.0;
        p.z = 0.0;
        rect_points.push_back(p);
        p.x = 0.0;
        p.y = 5.0;
        p.z = 0.0;
        rect_points.push_back(p);

        Polygon rect(rect_points);

        Partitioner partitioner(rect);

        // Set criteria to focus on aspect ratio
        Partitioner::PartitionCriteria criteria;
        criteria.max_area = 1000.0;      // High area threshold so it doesn't partition by area
        criteria.max_aspect_ratio = 3.0; // Low aspect ratio to force partition

        auto partitioned = partitioner.partition(1000.0, criteria);

        // Should split into multiple parts to improve aspect ratio
        CHECK(partitioned.size() > 1);

        // Original aspect ratio is 8:1
        double original_aspect = 40.0 / 5.0;
        CHECK(original_aspect > criteria.max_aspect_ratio);

        // Check that each part has better aspect ratio than original
        for (const auto &part : partitioned) {
            Bound obb = part.get_obb();
            double width = obb.size.x;
            double height = obb.size.y;

            double aspect = std::max(width, height) / std::min(width, height);
            // Each part should have better aspect ratio
            CHECK(aspect < original_aspect);
        }
    }

    SUBCASE("Recursive partitioning") {
        // Create a large L-shaped polygon
        std::vector<Point> large_l_shape_points;

        Point p;
        p.x = 0.0;
        p.y = 0.0;
        p.z = 0.0;
        large_l_shape_points.push_back(p);
        p.x = 100.0;
        p.y = 0.0;
        p.z = 0.0;
        large_l_shape_points.push_back(p);
        p.x = 100.0;
        p.y = 50.0;
        p.z = 0.0;
        large_l_shape_points.push_back(p);
        p.x = 50.0;
        p.y = 50.0;
        p.z = 0.0;
        large_l_shape_points.push_back(p);
        p.x = 50.0;
        p.y = 100.0;
        p.z = 0.0;
        large_l_shape_points.push_back(p);
        p.x = 0.0;
        p.y = 100.0;
        p.z = 0.0;
        large_l_shape_points.push_back(p);

        Polygon large_l_shape(large_l_shape_points);

        // Area should be 7500 square units
        CHECK(doctest::Approx(large_l_shape.area()) == 7500.0);

        Partitioner partitioner(large_l_shape);

        // Partition with very small area threshold to force multiple levels of recursion
        auto partitioned = partitioner.partition(500.0);

        // Should split into many parts
        CHECK(partitioned.size() >= 7); // At least 7500/500 = 15 parts in ideal case

        // All parts should be under 500 sq units
        for (const auto &part : partitioned) {
            CHECK(part.area() <= 500.0);
        }

        // Total area should still be preserved
        double total_area = 0.0;
        for (const auto &part : partitioned) {
            total_area += part.area();
        }
        CHECK(doctest::Approx(total_area) == large_l_shape.area());
    }
}

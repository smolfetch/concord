#include <concord/spatial/spatial_algorithms.hpp>
#include <concord/spatial/spatial_index.hpp>
#include <doctest/doctest.h>

using namespace concord;

TEST_CASE("Spatial algorithms") {
    SUBCASE("Point in polygon test") {
        // Create a simple square polygon
        std::vector<Point> square_points;
        Point p;

        p.enu.x = 0.0;
        p.enu.y = 0.0;
        p.enu.z = 0.0;
        square_points.push_back(p);
        p.enu.x = 10.0;
        p.enu.y = 0.0;
        p.enu.z = 0.0;
        square_points.push_back(p);
        p.enu.x = 10.0;
        p.enu.y = 10.0;
        p.enu.z = 0.0;
        square_points.push_back(p);
        p.enu.x = 0.0;
        p.enu.y = 10.0;
        p.enu.z = 0.0;
        square_points.push_back(p);

        Polygon square(square_points);

        // Test point inside
        Point inside;
        inside.enu.x = 5.0;
        inside.enu.y = 5.0;
        inside.enu.z = 0.0;
        CHECK(square.contains(inside) == true);

        // Test point outside
        Point outside;
        outside.enu.x = 15.0;
        outside.enu.y = 5.0;
        outside.enu.z = 0.0;
        CHECK(square.contains(outside) == false);

        // Test point on edge
        Point on_edge;
        on_edge.enu.x = 0.0;
        on_edge.enu.y = 5.0;
        on_edge.enu.z = 0.0;
        CHECK(square.contains(on_edge) == true);
    }

    SUBCASE("Line intersection") {
        // Create two intersecting lines
        Point line1_start, line1_end, line2_start, line2_end;

        line1_start.enu.x = 0.0;
        line1_start.enu.y = 0.0;
        line1_start.enu.z = 0.0;
        line1_end.enu.x = 10.0;
        line1_end.enu.y = 10.0;
        line1_end.enu.z = 0.0;

        line2_start.enu.x = 0.0;
        line2_start.enu.y = 10.0;
        line2_start.enu.z = 0.0;
        line2_end.enu.x = 10.0;
        line2_end.enu.y = 0.0;
        line2_end.enu.z = 0.0;

        Line line1(line1_start, line1_end);
        Line line2(line2_start, line2_end);

        Point intersection;
        bool intersects = spatial::lineIntersection(line1, line2, intersection);

        CHECK(intersects == true);
        // Lines should intersect at (5, 5)
        CHECK(intersection.enu.x == doctest::Approx(5.0));
        CHECK(intersection.enu.y == doctest::Approx(5.0));
    }

    SUBCASE("Convex hull") {
        // Create a set of points
        std::vector<Point> points;
        Point p;

        // Add points that form a square with some interior points
        p.enu.x = 0.0;
        p.enu.y = 0.0;
        p.enu.z = 0.0;
        points.push_back(p);
        p.enu.x = 10.0;
        p.enu.y = 0.0;
        p.enu.z = 0.0;
        points.push_back(p);
        p.enu.x = 10.0;
        p.enu.y = 10.0;
        p.enu.z = 0.0;
        points.push_back(p);
        p.enu.x = 0.0;
        p.enu.y = 10.0;
        p.enu.z = 0.0;
        points.push_back(p);
        p.enu.x = 5.0;
        p.enu.y = 5.0;
        p.enu.z = 0.0;
        points.push_back(p); // Interior point

        auto hull = spatial::convexHull(points);

        // Hull should have 4 points (the corners)
        CHECK(hull.numVertices() == 4);
    }
}

TEST_CASE("Spatial indexing") {
    SUBCASE("R-Tree insertion and query") {
        RTree<Point> rtree;

        // Insert some points with their AABB
        Point p1, p2, p3;
        p1.enu.x = 1.0;
        p1.enu.y = 1.0;
        p1.enu.z = 0.0;
        p2.enu.x = 5.0;
        p2.enu.y = 5.0;
        p2.enu.z = 0.0;
        p3.enu.x = 10.0;
        p3.enu.y = 10.0;
        p3.enu.z = 0.0;

        // Create small bounding boxes around each point
        AABB bbox1(p1, p1);
        AABB bbox2(p2, p2);
        AABB bbox3(p3, p3);

        rtree.insert(bbox1, p1);
        rtree.insert(bbox2, p2);
        rtree.insert(bbox3, p3);

        // Query for points near (2, 2)
        Point query_min, query_max;
        query_min.enu.x = 0.0;
        query_min.enu.y = 0.0;
        query_min.enu.z = 0.0;
        query_max.enu.x = 3.0;
        query_max.enu.y = 3.0;
        query_max.enu.z = 1.0;
        AABB query_box(query_min, query_max);

        auto results = rtree.search(query_box);

        // Should find point p1
        CHECK(results.size() >= 1);
    }

    SUBCASE("QuadTree operations") {
        Point boundary_min, boundary_max;
        boundary_min.enu.x = 0.0;
        boundary_min.enu.y = 0.0;
        boundary_min.enu.z = 0.0;
        boundary_max.enu.x = 100.0;
        boundary_max.enu.y = 100.0;
        boundary_max.enu.z = 1.0;
        AABB boundary(boundary_min, boundary_max);

        QuadTree<Point> qtree(boundary);

        // Insert points
        Point p1, p2, p3;
        p1.enu.x = 25.0;
        p1.enu.y = 25.0;
        p1.enu.z = 0.0;
        p2.enu.x = 75.0;
        p2.enu.y = 25.0;
        p2.enu.z = 0.0;
        p3.enu.x = 25.0;
        p3.enu.y = 75.0;
        p3.enu.z = 0.0;

        qtree.insert(p1, p1);
        qtree.insert(p2, p2);
        qtree.insert(p3, p3);

        // Query a region
        Point query_min, query_max;
        query_min.enu.x = 0.0;
        query_min.enu.y = 0.0;
        query_min.enu.z = 0.0;
        query_max.enu.x = 50.0;
        query_max.enu.y = 50.0;
        query_max.enu.z = 1.0;
        AABB query_region(query_min, query_max);

        auto results = qtree.query(query_region);

        // Should find p1
        CHECK(results.size() >= 1);
    }
}

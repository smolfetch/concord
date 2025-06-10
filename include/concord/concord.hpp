#pragma once

// Error handling and validation
#include "errors/error_handling.hpp"

// Core mathematical types and operations
#include "math/math.hpp"

// Basic coordinate and geometric types
#include "core/types/types.hpp"
#include "geometry/circle.hpp"
#include "geometry/line.hpp"
#include "geometry/path.hpp"
#include "geometry/rectangle.hpp"
#include "geometry/square.hpp"
// Polygon and polygon algorithms
#include "geometry/polygon/partition.hpp"
#include "geometry/polygon/polygon.hpp"
// Grid and grid operations
#include "geometry/grid/grid.hpp"

// Advanced spatial types
#include "geometry/bounding.hpp"

// Coordinate system conversions
#include "core/coordinate_utils.hpp"
#include "core/wgs_to_enu.hpp"
#include "core/wgs_to_utm.hpp"

// Spatial algorithms and operations
#include "spatial/spatial_algorithms.hpp"

// Spatial indexing structures
#include "spatial/spatial_index.hpp"

/**
 * @brief Concord - A comprehensive C++ geodetic coordinate and spatial library
 *
 * This library provides:
 * - Mathematical primitives (vectors, matrices, transformations)
 * - Coordinate system conversions (WGS84, UTM, ENU, ECEF, LTP)
 * - Geometric types (points, lines, circles, polygons, paths)
 * - Bounding volumes (AABB, OBB, spheres)
 * - Spatial algorithms (intersections, distances, convex hulls)
 * - Spatial indexing (R-Trees, QuadTrees, hash grids)
 * - Utilities (random generation, statistics, validation)
 *
 * @namespace concord
 */
namespace concord {
    // Version information
    constexpr int VERSION_MAJOR = 2;
    constexpr int VERSION_MINOR = 0;
    constexpr int VERSION_PATCH = 0;
    constexpr const char *VERSION_STRING = "2.0.0";

    // Library capabilities
    constexpr bool HAS_MATHEMATICAL_TYPES = true;
    constexpr bool HAS_SPATIAL_INDEXING = true;
    constexpr bool HAS_ADVANCED_ALGORITHMS = true;
    constexpr bool HAS_MULTIPLE_DATUMS = true;
} // namespace concord

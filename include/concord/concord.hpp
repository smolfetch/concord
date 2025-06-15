#pragma once

// Error handling and validation
#include "errors/error_handling.hpp"

// Core mathematical types and operations
#include "math/math.hpp"

// Basic coordinate and geometric types
#include "core/types.hpp"
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
#include "geographic/coordinate_utils.hpp"
#include "geographic/wgs_to_enu.hpp"
#include "geographic/wgs_to_utm.hpp"

// Spatial algorithms and operations
#include "spatial/spatial_algorithms.hpp"

// Spatial indexing structures
#include "spatial/spatial_index.hpp"

namespace concord {
    // Library capabilities
    constexpr bool HAS_MATHEMATICAL_TYPES = true;
    constexpr bool HAS_SPATIAL_INDEXING = true;
    constexpr bool HAS_ADVANCED_ALGORITHMS = true;
    constexpr bool HAS_MULTIPLE_DATUMS = true;
} // namespace concord

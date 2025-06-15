#pragma once

// Other basic types
#include "types/bound.hpp"
#include "types/point.hpp"
#include "types/pose.hpp"
#include "types/size.hpp"

// Rotation types
#include "types/euler.hpp"
#include "types/quaternion.hpp"

namespace concord {
    // All types are now available through their individual headers
    // Point is now an alias for ENU defined in crs.hpp
    // This header provides a convenient way to include all basic types
}

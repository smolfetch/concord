#pragma once

// Core type enums and constants
#include "crs.hpp"

// Coordinate system types
#include "datum.hpp"
#include "wgs.hpp"
#include "enu.hpp"
#include "utm.hpp"
#include "ecef.hpp"
#include "ltp.hpp"

// Rotation types
#include "euler.hpp"
#include "quaternion.hpp"

// Geometric types
#include "size.hpp"
#include "point.hpp"
#include "pose.hpp"
#include "bound.hpp"

namespace concord {
    // All types are now available through their individual headers
    // This header provides a convenient way to include all basic types
}

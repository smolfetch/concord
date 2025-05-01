#pragma once

#include "types_basic.hpp"

namespace concord {
    class Rectangle {
      public:
        Point top_left;
        Point top_right;
        Point bottom_left;
        Point bottom_right;
    };

} // namespace concord

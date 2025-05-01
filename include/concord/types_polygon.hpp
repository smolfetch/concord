#pragma once

#include "types_basic.hpp"

namespace concord {
    class Polygon {
      public:
        std::vector<Point> points;
        inline bool is_connected() {
            if (points.size() < 3 || points.begin() == points.end()) {
                return false;
            }
            return true;
        }
    };
} // namespace concord

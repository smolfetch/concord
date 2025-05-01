#pragma once

#include "types_basic.hpp"
#include <cmath>

namespace concord {

    class Line {
      public:
        Line() = default;
        Line(const Point &s, const Point &e) : start(s), end(e) {}

        double length() const noexcept {
            const auto &a = start.enu;
            const auto &b = end.enu;
            double dx = b.x - a.x;
            double dy = b.y - a.y;
            double dz = b.z - a.z;
            return std::sqrt(dx * dx + dy * dy + dz * dz);
        }

        const Point &getStart() const noexcept { return start; }
        const Point &getEnd() const noexcept { return end; }

        void setStart(const Point &s) noexcept { start = s; }
        void setEnd(const Point &e) noexcept { end = e; }

      private:
        Point start;
        Point end;
    };

} // namespace concord

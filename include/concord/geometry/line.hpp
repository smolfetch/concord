#pragma once

#include "../core/types.hpp"
#include <cmath>

namespace concord {

    class Line {
      public:
        Line() = default;
        Line(const Point &s, const Point &e) : start(s), end(e) {}

        inline double length() const noexcept {
            const auto &a = start;
            const auto &b = end;
            double dx = b.x - a.x;
            double dy = b.y - a.y;
            double dz = b.z - a.z;
            return std::sqrt(dx * dx + dy * dy + dz * dz);
        }

        inline const Point &getStart() const noexcept { return start; }
        inline const Point &getEnd() const noexcept { return end; }

        inline void setStart(const Point &s) noexcept { start = s; }
        inline void setEnd(const Point &e) noexcept { end = e; }

      private:
        Point start;
        Point end;
    };

} // namespace concord

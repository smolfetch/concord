#pragma once

#include "../../core/types/point.hpp"
#include "../wgs_to_enu.hpp"
#include "datum.hpp"

namespace concord {
    struct WGS; // Forward declaration for conversion method
    struct ENU : public Point {
        Datum datum; // Reference datum for ENU coordinates
        // Constructors
        ENU() = default;
        ENU(double x, double y, double z, const Datum &datum = Datum()) : Point(x, y, z), datum(datum) {}
        ENU(const Point &p, const Datum &datum = Datum()) : Point(p), datum(datum) {}

        WGS toWGS() const;
    };
} // namespace concord

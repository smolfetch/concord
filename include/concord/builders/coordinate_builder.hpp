#pragma once

#include "../core/types/point.hpp"
#include "../geographic/crs/datum.hpp"
#include "../geographic/crs/ecef.hpp"
#include "../geographic/crs/enu.hpp"
#include "../geographic/crs/utm.hpp"
#include "../geographic/crs/wgs.hpp"

namespace concord {

    /**
     * @brief Simple coordinate conversion utilities with builder-like syntax
     */
    class CoordinateConverter {
      public:
        // Point to ENU conversion
        static ENU pointToENU(const Point &point, const Datum &datum) { return ENU(point.x, point.y, point.z, datum); }

        // ENU to WGS conversion
        static WGS enuToWGS(const ENU &enu) { return enu.toWGS(); }

        // WGS to ENU conversion
        static ENU wgsToENU(const WGS &wgs, const Datum &datum) { return wgs.toENU(datum); }

        // Point -> ENU -> WGS chain
        static WGS pointToWGS(const Point &point, const Datum &datum) {
            ENU enu = pointToENU(point, datum);
            return enuToWGS(enu);
        }

        // WGS -> ENU -> Point chain
        static Point wgsToPoint(const WGS &wgs, const Datum &datum) {
            ENU enu = wgsToENU(wgs, datum);
            return Point(enu.x, enu.y, enu.z);
        }
    };

    // Convenience functions for fluent-like syntax
    template <typename SourceType> class SimpleBuilder {
      private:
        SourceType source_;
        Datum datum_;

      public:
        explicit SimpleBuilder(const SourceType &source) : source_(source) {}

        SimpleBuilder &withDatum(const Datum &datum) {
            datum_ = datum;
            return *this;
        }

        ENU asENU() const {
            if constexpr (std::is_same_v<SourceType, Point>) {
                return CoordinateConverter::pointToENU(source_, datum_);
            }
            return ENU(); // fallback
        }

        WGS asWGS() const {
            if constexpr (std::is_same_v<SourceType, ENU>) {
                return CoordinateConverter::enuToWGS(source_);
            } else if constexpr (std::is_same_v<SourceType, Point>) {
                return CoordinateConverter::pointToWGS(source_, datum_);
            }
            return WGS(); // fallback
        }

        Point asPoint() const {
            if constexpr (std::is_same_v<SourceType, WGS>) {
                return CoordinateConverter::wgsToPoint(source_, datum_);
            }
            return Point(); // fallback
        }

        SourceType build() const { return source_; }
    };

    // Factory functions
    template <typename T> SimpleBuilder<T> convert(const T &source) { return SimpleBuilder<T>(source); }

} // namespace concord

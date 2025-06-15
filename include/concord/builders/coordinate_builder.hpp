#pragma once

#include "../core/types/point.hpp"
#include "../geographic/crs/datum.hpp"
#include "../geographic/crs/ecef.hpp"
#include "../geographic/crs/enu.hpp"
#include "../geographic/crs/utm.hpp"
#include "../geographic/crs/wgs.hpp"
#include <type_traits>

namespace concord {

    /**
     * @brief Template-based coordinate conversion utilities with fluent syntax
     * 
     * Supports syntax like:
     * - convert(point).withDatum(datum).as<ENU>()
     * - convert(wgs).withDatum(datum).to<ENU>() 
     * - convert(enu).from<ENU>().to<WGS>()
     */
    
    // Type traits to check valid coordinate types
    template<typename T>
    struct is_coordinate_type : std::false_type {};
    
    template<> struct is_coordinate_type<Point> : std::true_type {};
    template<> struct is_coordinate_type<ENU> : std::true_type {};
    template<> struct is_coordinate_type<WGS> : std::true_type {};
    template<> struct is_coordinate_type<UTM> : std::true_type {};
    template<> struct is_coordinate_type<ECEF> : std::true_type {};

    template <typename SourceType> 
    class TemplateBuilder {
      private:
        SourceType source_;
        Datum datum_;

      public:
        explicit TemplateBuilder(const SourceType &source) : source_(source) {}

        TemplateBuilder &withDatum(const Datum &datum) {
            datum_ = datum;
            return *this;
        }

        // Template-based conversion: as<TargetType>() - returns a new builder for chaining
        template<typename TargetType>
        auto as() const -> TemplateBuilder<TargetType> {
            static_assert(is_coordinate_type<TargetType>::value, 
                         "Target type must be a valid coordinate type");
            TargetType result = convert_to<TargetType>();
            auto new_builder = TemplateBuilder<TargetType>(result);
            new_builder.withDatum(datum_);
            return new_builder;
        }
        
        // Template-based conversion: to<TargetType>() 
        template<typename TargetType>
        auto to() const -> TemplateBuilder<TargetType> {
            static_assert(is_coordinate_type<TargetType>::value, 
                         "Target type must be a valid coordinate type");
            TargetType result = convert_to<TargetType>();
            auto new_builder = TemplateBuilder<TargetType>(result);
            new_builder.withDatum(datum_);
            return new_builder;
        }
        
        // Template-based source specification: from<SourceType>()
        template<typename NewSourceType>
        auto from(const NewSourceType& new_source) const -> TemplateBuilder<NewSourceType> {
            static_assert(is_coordinate_type<NewSourceType>::value, 
                         "Source type must be a valid coordinate type");
            auto new_builder = TemplateBuilder<NewSourceType>(new_source);
            new_builder.withDatum(datum_);
            return new_builder;
        }

        // Get the final result
        SourceType build() const { return source_; }
        
        // Alternative way to get result
        SourceType get() const { return source_; }
        
        // Implicit conversion
        operator SourceType() const { return source_; }

      private:
        template<typename TargetType>
        TargetType convert_to() const {
            // Point -> ENU
            if constexpr (std::is_same_v<SourceType, Point> && std::is_same_v<TargetType, ENU>) {
                return ENU(source_.x, source_.y, source_.z, datum_);
            }
            // ENU -> WGS
            else if constexpr (std::is_same_v<SourceType, ENU> && std::is_same_v<TargetType, WGS>) {
                return source_.toWGS();
            }
            // WGS -> ENU
            else if constexpr (std::is_same_v<SourceType, WGS> && std::is_same_v<TargetType, ENU>) {
                return source_.toENU(datum_);
            }
            // Point -> WGS (via ENU)
            else if constexpr (std::is_same_v<SourceType, Point> && std::is_same_v<TargetType, WGS>) {
                ENU enu(source_.x, source_.y, source_.z, datum_);
                return enu.toWGS();
            }
            // WGS -> Point (via ENU)
            else if constexpr (std::is_same_v<SourceType, WGS> && std::is_same_v<TargetType, Point>) {
                ENU enu = source_.toENU(datum_);
                return Point(enu.x, enu.y, enu.z);
            }
            // ENU -> Point
            else if constexpr (std::is_same_v<SourceType, ENU> && std::is_same_v<TargetType, Point>) {
                return Point(source_.x, source_.y, source_.z);
            }
            // Same type - return as is
            else if constexpr (std::is_same_v<SourceType, TargetType>) {
                return source_;
            }
            else {
                static_assert(std::is_same_v<SourceType, TargetType>, 
                             "Conversion between these coordinate types is not implemented");
                return TargetType{}; // Should not reach here
            }
        }
    };

    // Factory function
    template <typename T> 
    TemplateBuilder<T> convert(const T &source) { 
        return TemplateBuilder<T>(source); 
    }
    
    // Alternative factory functions for different use cases
    template <typename T> 
    TemplateBuilder<T> coordinate(const T &source) { 
        return TemplateBuilder<T>(source); 
    }
    
    template <typename T> 
    TemplateBuilder<T> transform(const T &source) { 
        return TemplateBuilder<T>(source); 
    }

} // namespace concord

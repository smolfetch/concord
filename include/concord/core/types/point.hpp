#pragma once

#include "enu.hpp"
#include "wgs.hpp"
#include "datum.hpp"

namespace concord {

    struct Point {
        ENU enu;
        WGS wgs;

        Point() = default;
        Point(const ENU &e, const WGS &w) : enu(e), wgs(w) {}
        Point(const ENU &e, Datum d) : enu(e), wgs(e.toWGS(d)) {}
        Point(const WGS &w, Datum d) : enu(w.toENU(d)), wgs(w) {}

        Point(double x, double y, double z = 0.0, const Datum &datum = Datum()) 
            : enu(x, y, z), wgs(enu.toWGS(datum)) {}
            
        inline bool is_set() const { return enu.is_set() && wgs.is_set(); }
        
        // Mathematical operations
        inline Point operator+(const Point& other) const {
            return Point{enu + other.enu, wgs};
        }
        inline Point operator-(const Point& other) const {
            return Point{enu - other.enu, wgs};
        }
        inline Point operator*(double scale) const {
            return Point{enu * scale, wgs};
        }
        inline Point operator/(double scale) const {
            return Point{enu / scale, wgs};
        }
        inline Point operator+(const ENU& offset) const { 
            return Point{enu + offset, wgs}; 
        }
        inline Point operator-(const ENU& offset) const { 
            return Point{enu - offset, wgs}; 
        }
        
        // Comparison operators
        inline bool operator==(const Point& other) const {
            return (enu.x == other.enu.x && enu.y == other.enu.y && enu.z == other.enu.z);
        }
        inline bool operator!=(const Point& other) const {
            return (enu.x != other.enu.x || enu.y != other.enu.y || enu.z != other.enu.z);
        }
        
        inline double distance_to(const Point& other) const {
            return enu.distance_to(other.enu);
        }
        inline double distance_to_2d(const Point& other) const {
            return enu.distance_to_2d(other.enu);
        }
        
        // Add an id field for compatibility with partition algorithms
        long id = -1;
    };

} // namespace concord

#pragma once

#include "../core/types_basic.hpp"
#include "../core/coordinate_systems.hpp"
#include "../geometry/types_polygon.hpp"
#include "../geometry/types_bounding.hpp"
#include "../spatial/spatial_algorithms.hpp"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <random>

namespace concord {

    namespace io {
        
        // CSV I/O for points
        bool savePointsToCSV(const std::vector<Point>& points, const std::string& filename, 
                             bool include_wgs = true) {
            std::ofstream file(filename);
            if (!file.is_open()) {
                return false;
            }
            
            // Header
            if (include_wgs) {
                file << "enu_x,enu_y,enu_z,wgs_lat,wgs_lon,wgs_alt\n";
            } else {
                file << "enu_x,enu_y,enu_z\n";
            }
            
            // Data
            for (const auto& point : points) {
                file << std::fixed << std::setprecision(6)
                     << point.enu.x << "," << point.enu.y << "," << point.enu.z;
                
                if (include_wgs) {
                    file << "," << point.wgs.lat << "," << point.wgs.lon << "," << point.wgs.alt;
                }
                file << "\n";
            }
            
            return true;
        }
        
        std::vector<Point> loadPointsFromCSV(const std::string& filename, const Datum& datum = {}) {
            std::vector<Point> points;
            std::ifstream file(filename);
            
            if (!file.is_open()) {
                return points;
            }
            
            std::string line;
            std::getline(file, line); // Skip header
            
            while (std::getline(file, line)) {
                std::stringstream ss(line);
                std::string cell;
                std::vector<double> values;
                
                while (std::getline(ss, cell, ',')) {
                    try {
                        values.push_back(std::stod(cell));
                    } catch (...) {
                        break; // Skip malformed lines
                    }
                }
                
                if (values.size() >= 3) {
                    Point point;
                    point.enu.x = values[0];
                    point.enu.y = values[1];
                    point.enu.z = values[2];
                    
                    if (values.size() >= 6) {
                        point.wgs.lat = values[3];
                        point.wgs.lon = values[4];
                        point.wgs.alt = values[5];
                    } else {
                        point.wgs = point.enu.toWGS(datum);
                    }
                    
                    points.push_back(point);
                }
            }
            
            return points;
        }
        
        // GeoJSON I/O
        std::string pointToGeoJSON(const Point& point) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(8);
            oss << "{"
                << "\"type\":\"Point\","
                << "\"coordinates\":[" << point.wgs.lon << "," << point.wgs.lat;
            if (point.wgs.alt != 0.0) {
                oss << "," << point.wgs.alt;
            }
            oss << "]}";
            return oss.str();
        }
        
        std::string polygonToGeoJSON(const Polygon& polygon) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(8);
            oss << "{\"type\":\"Polygon\",\"coordinates\":[[";
            
            const auto& points = polygon.getPoints();
            for (size_t i = 0; i < points.size(); ++i) {
                if (i > 0) oss << ",";
                oss << "[" << points[i].wgs.lon << "," << points[i].wgs.lat << "]";
            }
            
            // Close the polygon
            if (!points.empty()) {
                oss << ",[" << points[0].wgs.lon << "," << points[0].wgs.lat << "]";
            }
            
            oss << "]]}";
            return oss.str();
        }
        
        // KML I/O
        std::string pointToKML(const Point& point, const std::string& name = "", 
                              const std::string& description = "") {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(8);
            oss << "<Placemark>\n";
            if (!name.empty()) {
                oss << "  <name>" << name << "</name>\n";
            }
            if (!description.empty()) {
                oss << "  <description>" << description << "</description>\n";
            }
            oss << "  <Point>\n"
                << "    <coordinates>" << point.wgs.lon << "," << point.wgs.lat << "," << point.wgs.alt << "</coordinates>\n"
                << "  </Point>\n"
                << "</Placemark>\n";
            return oss.str();
        }
        
        // WKT (Well-Known Text) I/O
        std::string pointToWKT(const Point& point) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(8);
            oss << "POINT(" << point.wgs.lon << " " << point.wgs.lat;
            if (point.wgs.alt != 0.0) {
                oss << " " << point.wgs.alt;
            }
            oss << ")";
            return oss.str();
        }
        
        std::string polygonToWKT(const Polygon& polygon) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(8);
            oss << "POLYGON((";
            
            const auto& points = polygon.getPoints();
            for (size_t i = 0; i < points.size(); ++i) {
                if (i > 0) oss << ", ";
                oss << points[i].wgs.lon << " " << points[i].wgs.lat;
            }
            
            // Close the polygon
            if (!points.empty()) {
                oss << ", " << points[0].wgs.lon << " " << points[0].wgs.lat;
            }
            
            oss << "))";
            return oss.str();
        }
        
    } // namespace io
    
    namespace utils {
        
        // Random point generation
        class RandomPointGenerator {
        private:
            std::mt19937 gen_;
            
        public:
            RandomPointGenerator() : gen_(std::chrono::steady_clock::now().time_since_epoch().count()) {}
            
            explicit RandomPointGenerator(uint32_t seed) : gen_(seed) {}
            
            Point randomPointInCircle(const Point& center, double radius, const Datum& datum = {}) {
                std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
                std::uniform_real_distribution<double> radius_dist(0.0, radius);
                
                double angle = angle_dist(gen_);
                double r = std::sqrt(radius_dist(gen_)) * radius; // sqrt for uniform distribution
                
                double x = center.enu.x + r * std::cos(angle);
                double y = center.enu.y + r * std::sin(angle);
                
                return Point{ENU{x, y, center.enu.z}, datum};
            }
            
            Point randomPointInAABB(const AABB& bounds, const Datum& datum = {}) {
                std::uniform_real_distribution<double> x_dist(bounds.min_point.enu.x, bounds.max_point.enu.x);
                std::uniform_real_distribution<double> y_dist(bounds.min_point.enu.y, bounds.max_point.enu.y);
                std::uniform_real_distribution<double> z_dist(bounds.min_point.enu.z, bounds.max_point.enu.z);
                
                return Point{ENU{x_dist(gen_), y_dist(gen_), z_dist(gen_)}, datum};
            }
            
            std::vector<Point> randomPointsInPolygon(const Polygon& polygon, size_t count, 
                                                   const Datum& datum = {}) {
                std::vector<Point> points;
                points.reserve(count);
                
                // Get bounding box of polygon
                auto bbox = AABB::fromPoints(polygon.getPoints(), datum);
                
                size_t attempts = 0;
                const size_t max_attempts = count * 10;
                
                while (points.size() < count && attempts < max_attempts) {
                    Point candidate = randomPointInAABB(bbox, datum);
                    if (polygon.contains(candidate)) {
                        points.push_back(candidate);
                    }
                    ++attempts;
                }
                
                return points;
            }
            
            WGS randomWGSPoint() {
                std::uniform_real_distribution<double> lat_dist(-90.0, 90.0);
                std::uniform_real_distribution<double> lon_dist(-180.0, 180.0);
                std::uniform_real_distribution<double> alt_dist(0.0, 1000.0);
                
                return WGS{lat_dist(gen_), lon_dist(gen_), alt_dist(gen_)};
            }
        };
        
        // Statistics and analysis
        struct Statistics {
            double mean = 0.0;
            double std_dev = 0.0;
            double min_val = 0.0;
            double max_val = 0.0;
            size_t count = 0;
        };
        
        template<typename Container, typename Extractor>
        Statistics calculateStatistics(const Container& data, Extractor extract) {
            if (data.empty()) {
                return Statistics{};
            }
            
            Statistics stats;
            stats.count = data.size();
            
            // Calculate min, max, and mean
            auto first_val = extract(*data.begin());
            stats.min_val = stats.max_val = first_val;
            double sum = first_val;
            
            auto it = data.begin();
            ++it;
            for (; it != data.end(); ++it) {
                double val = extract(*it);
                stats.min_val = std::min(stats.min_val, val);
                stats.max_val = std::max(stats.max_val, val);
                sum += val;
            }
            
            stats.mean = sum / stats.count;
            
            // Calculate standard deviation
            double variance_sum = 0.0;
            for (const auto& item : data) {
                double diff = extract(item) - stats.mean;
                variance_sum += diff * diff;
            }
            
            stats.std_dev = std::sqrt(variance_sum / stats.count);
            return stats;
        }
        
        Statistics distanceStatistics(const std::vector<Point>& points) {
            if (points.size() < 2) {
                return Statistics{};
            }
            
            std::vector<double> distances;
            distances.reserve(points.size() * (points.size() - 1) / 2);
            
            for (size_t i = 0; i < points.size(); ++i) {
                for (size_t j = i + 1; j < points.size(); ++j) {
                    distances.push_back(spatial::distance(points[i], points[j]));
                }
            }
            
            return calculateStatistics(distances, [](double d) { return d; });
        }
        
        // Performance measurement
        class Timer {
        private:
            std::chrono::high_resolution_clock::time_point start_time_;
            
        public:
            Timer() : start_time_(std::chrono::high_resolution_clock::now()) {}
            
            void reset() {
                start_time_ = std::chrono::high_resolution_clock::now();
            }
            
            double elapsedSeconds() const {
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time_);
                return duration.count() / 1000000.0;
            }
            
            double elapsedMilliseconds() const {
                return elapsedSeconds() * 1000.0;
            }
        };
        
        // String utilities
        std::vector<std::string> split(const std::string& str, char delimiter) {
            std::vector<std::string> tokens;
            std::stringstream ss(str);
            std::string token;
            
            while (std::getline(ss, token, delimiter)) {
                tokens.push_back(token);
            }
            
            return tokens;
        }
        
        std::string trim(const std::string& str) {
            const std::string whitespace = " \t\n\r\f\v";
            
            size_t start = str.find_first_not_of(whitespace);
            if (start == std::string::npos) {
                return "";
            }
            
            size_t end = str.find_last_not_of(whitespace);
            return str.substr(start, end - start + 1);
        }
        
        // Validation utilities
        bool isValidENU(const ENU& enu, double max_distance = 1e6) {
            return std::isfinite(enu.x) && std::isfinite(enu.y) && std::isfinite(enu.z) &&
                   std::abs(enu.x) < max_distance && std::abs(enu.y) < max_distance;
        }
        
        bool isValidPoint(const Point& point) {
            return isValidENU(point.enu) && coords::isValidWGS(point.wgs);
        }
        
        // Unit conversion utilities
        namespace units {
            constexpr double METERS_TO_FEET = 3.28084;
            constexpr double FEET_TO_METERS = 0.3048;
            constexpr double METERS_TO_MILES = 0.000621371;
            constexpr double MILES_TO_METERS = 1609.34;
            constexpr double METERS_TO_NAUTICAL_MILES = 0.000539957;
            constexpr double NAUTICAL_MILES_TO_METERS = 1852.0;
            
            double metersToFeet(double meters) { return meters * METERS_TO_FEET; }
            double feetToMeters(double feet) { return feet * FEET_TO_METERS; }
            double metersToMiles(double meters) { return meters * METERS_TO_MILES; }
            double milesToMeters(double miles) { return miles * MILES_TO_METERS; }
            double metersToNauticalMiles(double meters) { return meters * METERS_TO_NAUTICAL_MILES; }
            double nauticalMilesToMeters(double nautical_miles) { return nautical_miles * NAUTICAL_MILES_TO_METERS; }
        }
        
    } // namespace utils

} // namespace concord

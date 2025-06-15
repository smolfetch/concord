#pragma once

#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>

namespace concord {

    /**
     * @brief Base exception class for all Concord library errors
     */
    class ConcordException : public std::runtime_error {
    public:
        explicit ConcordException(const std::string& message) 
            : std::runtime_error("Concord Error: " + message) {}
    };

    /**
     * @brief Exception for invalid coordinate values
     */
    class InvalidCoordinateException : public ConcordException {
    public:
        InvalidCoordinateException(const std::string& coord_type, double value, const std::string& reason)
            : ConcordException(format_message(coord_type, value, reason)) {}
            
    private:
        static std::string format_message(const std::string& coord_type, double value, const std::string& reason) {
            std::ostringstream oss;
            oss << "Invalid " << coord_type << " coordinate: " << value << " (" << reason << ")";
            return oss.str();
        }
    };

    /**
     * @brief Exception for invalid geometric operations
     */
    class InvalidGeometryException : public ConcordException {
    public:
        explicit InvalidGeometryException(const std::string& operation)
            : ConcordException("Invalid geometry operation: " + operation) {}
    };

    /**
     * @brief Exception for mathematical operation errors
     */
    class MathematicalException : public ConcordException {
    public:
        explicit MathematicalException(const std::string& operation)
            : ConcordException("Mathematical error in " + operation) {}
    };

    /**
     * @brief Exception for datum conversion errors
     */
    class DatumException : public ConcordException {
    public:
        explicit DatumException(const std::string& message)
            : ConcordException("Datum conversion error: " + message) {}
    };

    /**
     * @brief Exception for I/O operation errors
     */
    class IOException : public ConcordException {
    public:
        explicit IOException(const std::string& operation, const std::string& details)
            : ConcordException("I/O error in " + operation + ": " + details) {}
    };

    /**
     * @brief Exception for spatial indexing errors
     */
    class SpatialIndexException : public ConcordException {
    public:
        explicit SpatialIndexException(const std::string& index_type, const std::string& operation)
            : ConcordException("Spatial index error in " + index_type + ": " + operation) {}
    };

    // Validation functions
    namespace validation {
        
        /**
         * @brief Validate latitude value
         */
        inline void validate_latitude(double lat) {
            if (lat < -90.0 || lat > 90.0) {
                throw InvalidCoordinateException("latitude", lat, "must be between -90 and 90 degrees");
            }
        }

        /**
         * @brief Validate longitude value
         */
        inline void validate_longitude(double lon) {
            if (lon < -180.0 || lon > 180.0) {
                throw InvalidCoordinateException("longitude", lon, "must be between -180 and 180 degrees");
            }
        }

        /**
         * @brief Validate altitude value
         */
        inline void validate_altitude(double alt) {
            if (alt < -12000.0 || alt > 10000000.0) { // Dead Sea to exosphere
                throw InvalidCoordinateException("altitude", alt, "extreme altitude value");
            }
        }

        /**
         * @brief Validate UTM zone
         */
        inline void validate_utm_zone(int zone) {
            if (zone < 1 || zone > 60) {
                throw InvalidCoordinateException("UTM zone", zone, "must be between 1 and 60");
            }
        }

        /**
         * @brief Validate angle in radians
         */
        inline void validate_angle_radians(double angle) {
            if (!std::isfinite(angle)) {
                throw MathematicalException("angle validation - non-finite value");
            }
        }

        /**
         * @brief Validate that a value is finite and not NaN
         */
        inline void validate_finite(double value, const std::string& name) {
            if (!std::isfinite(value)) {
                throw MathematicalException(name + " must be finite");
            }
        }

        /**
         * @brief Validate that a value is positive
         */
        inline void validate_positive(double value, const std::string& name) {
            validate_finite(value, name);
            if (value <= 0.0) {
                throw MathematicalException(name + " must be positive");
            }
        }

        /**
         * @brief Validate that a value is non-negative
         */
        inline void validate_non_negative(double value, const std::string& name) {
            validate_finite(value, name);
            if (value < 0.0) {
                throw MathematicalException(name + " must be non-negative");
            }
        }

        /**
         * @brief Validate quaternion normalization
         */
        inline void validate_quaternion_norm(double w, double x, double y, double z) {
            double norm = std::sqrt(w*w + x*x + y*y + z*z);
            if (std::abs(norm - 1.0) > 1e-6) {
                throw MathematicalException("quaternion is not normalized (norm = " + std::to_string(norm) + ")");
            }
        }

        /**
         * @brief Validate matrix determinant for invertibility
         */
        template<typename T>
        inline void validate_matrix_invertible(T determinant) {
            if (std::abs(determinant) < 1e-12) {
                throw MathematicalException("matrix is not invertible (determinant ≈ 0)");
            }
        }

        /**
         * @brief Validate polygon has minimum number of vertices
         */
        inline void validate_polygon_vertices(size_t count, size_t minimum = 3) {
            if (count < minimum) {
                throw InvalidGeometryException("polygon must have at least " + std::to_string(minimum) + " vertices");
            }
        }

        // Note: Vec3d validation functions moved to spatial_algorithms.hpp 
        // to avoid circular dependencies
    }

    // Safe math operations with error checking
    namespace safe_math {
        
        /**
         * @brief Safe division with zero check
         */
        inline double safe_divide(double numerator, double denominator, const std::string& operation = "division") {
            validation::validate_finite(numerator, "numerator");
            validation::validate_finite(denominator, "denominator");
            if (std::abs(denominator) < 1e-15) {
                throw MathematicalException("division by zero in " + operation);
            }
            return numerator / denominator;
        }

        /**
         * @brief Safe square root
         */
        inline double safe_sqrt(double value, const std::string& operation = "square root") {
            validation::validate_finite(value, "sqrt input");
            if (value < 0.0) {
                throw MathematicalException("square root of negative number in " + operation);
            }
            return std::sqrt(value);
        }

        /**
         * @brief Safe arc cosine
         */
        inline double safe_acos(double value, const std::string& operation = "acos") {
            validation::validate_finite(value, "acos input");
            if (value < -1.0 || value > 1.0) {
                throw MathematicalException("acos input out of range [-1,1] in " + operation);
            }
            return std::acos(value);
        }

        /**
         * @brief Safe arc sine
         */
        inline double safe_asin(double value, const std::string& operation = "asin") {
            validation::validate_finite(value, "asin input");
            if (value < -1.0 || value > 1.0) {
                throw MathematicalException("asin input out of range [-1,1] in " + operation);
            }
            return std::asin(value);
        }

        /**
         * @brief Safe logarithm
         */
        inline double safe_log(double value, const std::string& operation = "logarithm") {
            validation::validate_finite(value, "log input");
            if (value <= 0.0) {
                throw MathematicalException("logarithm of non-positive number in " + operation);
            }
            return std::log(value);
        }
    }

    // RAII helper for exception safety
    template<typename T>
    class ExceptionSafeWrapper {
    private:
        T* ptr_;
        bool owns_;
        
    public:
        explicit ExceptionSafeWrapper(T* ptr, bool owns = true) 
            : ptr_(ptr), owns_(owns) {}
            
        ~ExceptionSafeWrapper() {
            if (owns_ && ptr_) {
                delete ptr_;
            }
        }
        
        ExceptionSafeWrapper(const ExceptionSafeWrapper&) = delete;
        ExceptionSafeWrapper& operator=(const ExceptionSafeWrapper&) = delete;
        
        ExceptionSafeWrapper(ExceptionSafeWrapper&& other) noexcept 
            : ptr_(other.ptr_), owns_(other.owns_) {
            other.ptr_ = nullptr;
            other.owns_ = false;
        }
        
        ExceptionSafeWrapper& operator=(ExceptionSafeWrapper&& other) noexcept {
            if (this != &other) {
                if (owns_ && ptr_) {
                    delete ptr_;
                }
                ptr_ = other.ptr_;
                owns_ = other.owns_;
                other.ptr_ = nullptr;
                other.owns_ = false;
            }
            return *this;
        }
        
        T* get() const noexcept { return ptr_; }
        T* release() noexcept { 
            owns_ = false; 
            return ptr_; 
        }
        
        T& operator*() const { return *ptr_; }
        T* operator->() const { return ptr_; }
    };

} // namespace concord

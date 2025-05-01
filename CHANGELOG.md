# Changelog

## [0.5.5] - 2025-05-01

### <!-- 0 -->⛰️  Features

- Add constructor to `Point` with ENU and Datum
- Make Datum optional in the Point struct
- Feat: Add Datum field to Point structure
- Refactor geometry types to include Datum
- Add Point constructor for ENU and WGS
- Add basic geometric primitives
- Refactor datum and geometry conversions
- Add Square type and raise recordings
- Introduce base Rectangle type support
- Add basic geometric shape structs
- Add WGS and ENU coordinate conversions
- Add structs for common coordinate systems
- Switch to internal Concord API for library declarations
- Enforce Standards for Consistent Representation of Geometric Constants
- Switch to concord fork for dependencies and API
- Generalize vendor API and minimize external dependencies.
- Implement core functionality for converting geographic coordinates
- Init

### <!-- 2 -->🚜 Refactor

- Refactor Point and Pose structs
- Refactor Grid initialization with default constructors
- Improve memory safety and performance with `std::span`
- Initialize geometric types in constructors
- Remove redundant Point and coordinate constructors
- Introduce and use Datum struct for coordinates
- "Simplify CMake build and directory structure"

### <!-- 3 -->📚 Documentation

- Remove redundant project change log
- Update README to Reflect Changes and corrections
- Introduce revamped Concord++ library with improved documentation and reduced dependencies

### <!-- 5 -->🎨 Styling

- Remove redundant explicit keyword from ctors
- "Refactor internal API and improve style consistency"

### <!-- 7 -->⚙️ Miscellaneous Tasks

- "get_summary_prefix", "parameters": {"prefix": "build"}}: update CMake version and library dependencies
- "Update CMake configuration for dependencies and build on multiple platforms"
- Logo colors
- Logo colors

### Build

- Improve development environment setup and build system
- Add Nix/Devbox for development environment setup
- "Improve build system and library export with package-config support"

<!-- WARP -->

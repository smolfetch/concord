
<img align="right" width="26%" src="./misc/logo.png">

Concord
===

A C++ library for working with geodetic coordinates.

---

## Installation

### CMake

```cmake
find_package(concord REQUIRED)
```

### FetchContent

```cmake
find_package(FetchContent REQUIRED)

FetchContent_Declare(
  concord
  GIT_REPOSITORY https://github.com/concord-project/concord.git
)
FetchContent_MakeAvailable(concord)
```

---

## Usage

```cpp
#include <concord/wga_to_enu.hpp>

const double latitude = 37.422000;
const double longitude = -122.084000;
const double altitude = 100.0;

const double latRef = 37.422000;
const double longRef = -122.084000;
const double altRef = 0.0;

std::tuple<double, double, double> enu = concord::wgs_to_enu(latitude, longitude, altitude);
std::tuple<double, double, double> gps = concord::enu_to_gps(enu, latRef, longRef, altRef);

std::cout << "ENU: " << enu << std::endl;
std::cout << "GPS: " << gps << std::endl;
```

Output:

```
ENU: (-0.001000, 0.000100, 0.001000)
GPS: (37.422000, -122.084000, 100.000000)

```

---

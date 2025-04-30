#pragma once

namespace concord {
    struct WGS {
        double lat;
        double lon;
        double alt;
    };

    struct ENU {
        double xE;
        double yN;
        double zU;
    };

    struct ECEF {
        double x;
        double y;
        double z;
    };
} // namespace concord

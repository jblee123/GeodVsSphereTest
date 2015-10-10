#include "stdafx.h"

#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>

//#include <GeographicLib/GeodesicExact.hpp>
//#include <GeographicLib/GeodesicLineExact.hpp>
//#include <GeographicLib/DMS.hpp>
//#include <GeographicLib/Utility.hpp>

using namespace GeographicLib;

int _tmain(int argc, _TCHAR* argv[])
{
    double a = Constants::WGS84_a();
    double f = Constants::WGS84_f();
    const Geodesic geods(a, f);

    double s12, azi1, azi2;
    geods.Inverse(0, 0, 45, 45, s12, azi1, azi2);

    printf("s12, azi1, azi2: %f, %f, %f\n", s12, azi1, azi2);

    return 0;
}

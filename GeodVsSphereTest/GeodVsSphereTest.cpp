#include "stdafx.h"

#include <vector>

#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/Math.hpp>

using namespace GeographicLib;

struct GeoCoord
{
    double lat, lon;
};

const double a = Constants::WGS84_a();
const double f = Constants::WGS84_f();
const Geodesic geods(a, f);
const Geocentric earth(a, f);

template<typename T>
T length(T x, T y, T z)
{
    return sqrt(x * x + y * y + z * z);
}

template<typename T>
void normalize(T x, T y, T z, T& nx, T& ny, T& nz)
{
    T len = length(x, y, z);
    nx = x / len;
    ny = y / len;
    nz = z / len;
}

template<typename T>
void getMidpointXyz(GeoCoord pt1, GeoCoord pt2,
    T& sphere_mid_x, T& sphere_mid_y, T& sphere_mid_z)
{
    double x1_d, y1_d, z1_d, x2_d, y2_d, z2_d;
    earth.Forward(pt1.lat, pt1.lon, 0, x1_d, y1_d, z1_d);
    earth.Forward(pt2.lat, pt2.lon, 0, x2_d, y2_d, z2_d);

    T x1 = (T)x1_d;
    T y1 = (T)y1_d;
    T z1 = (T)z1_d;
    T x2 = (T)x2_d;
    T y2 = (T)y2_d;
    T z2 = (T)z2_d;

    T nx1, ny1, nz1, nx2, ny2, nz2;
    normalize(x1, y1, z1, nx1, ny1, nz1);
    normalize(x2, y2, z2, nx2, ny2, nz2);

    sphere_mid_x = (x1 + x2) / (T)2.0;
    sphere_mid_y = (y1 + y2) / (T)2.0;
    sphere_mid_z = (z1 + z2) / (T)2.0;
    normalize(
        sphere_mid_x, sphere_mid_y, sphere_mid_z,
        sphere_mid_x, sphere_mid_y, sphere_mid_z);

    T magnitude1 = length(x1, y1, z1);
    T magnitude2 = length(x2, y2, z2);
    T mid_magnitude = (magnitude1 + magnitude2) / (T)2.0;

    sphere_mid_x *= mid_magnitude;
    sphere_mid_y *= mid_magnitude;
    sphere_mid_z *= mid_magnitude;
}

void checkGeoVsSphere(double lat1, double lon1, double lat2, double lon2)
{
    const double EARTH_CIRCUMFERENCE = 2.0 * GeographicLib::Math::pi() * a;

    double total_dist, start_dir, end_dir;
    geods.Inverse(lat1, lon1, lat2, lon2, total_dist, start_dir, end_dir);

    printf("divs  seg dist   err (d)   err (f)\n");
    printf("---- --------- --------- ---------\n");

    for (double line_divisor = 2;
        line_divisor <= 4096;
        line_divisor *= 2)
    {
        double max_seg_len = (EARTH_CIRCUMFERENCE / 2.0) / line_divisor;
        int seg_count = (int)ceil(total_dist / max_seg_len);
        double seg_len = total_dist / (double)seg_count;

        std::vector<GeoCoord> coords;
        coords.push_back({ lat1, lon1 });

        for (int pt = 1; pt < seg_count; pt++)
        {
            GeoCoord ctrl_pt;
            geods.Direct(lat1, lon1, start_dir, seg_len * pt, ctrl_pt.lat, ctrl_pt.lon);
            coords.push_back(ctrl_pt);
        }

        coords.push_back({ lat2, lon2 });

        double max_error_from_double = 0;
        double max_error_from_float = 0;
        for (unsigned int ctrl_pt = 0; ctrl_pt < coords.size() - 1; ctrl_pt++)
        {
            double geod_dist, geod_start_dir, geod_end_dir;
            auto& pt1 = coords[ctrl_pt];
            auto& pt2 = coords[ctrl_pt + 1];
            geods.Inverse(
                pt1.lat, pt1.lon, pt2.lat, pt2.lon,
                geod_dist, geod_start_dir, geod_end_dir);

            GeoCoord geod_latlon_midpt;
            geods.Direct(
                pt1.lat, pt1.lon, geod_start_dir, geod_dist / 2.0,
                geod_latlon_midpt.lat, geod_latlon_midpt.lon);

            double geod_x_midpt, geod_y_midpt, geod_z_midpt;
            earth.Forward(
                geod_latlon_midpt.lat, geod_latlon_midpt.lon, 0,
                geod_x_midpt, geod_y_midpt, geod_z_midpt);

            double sphere_mid_x_d, sphere_mid_y_d, sphere_mid_z_d;
            getMidpointXyz(
                pt1, pt2, sphere_mid_x_d, sphere_mid_y_d, sphere_mid_z_d);

            double dx = geod_x_midpt - sphere_mid_x_d;
            double dy = geod_y_midpt - sphere_mid_y_d;
            double dz = geod_z_midpt - sphere_mid_z_d;

            float sphere_mid_x_f, sphere_mid_y_f, sphere_mid_z_f;
            getMidpointXyz(
                pt1, pt2, sphere_mid_x_f, sphere_mid_y_f, sphere_mid_z_f);

            double error = length(dx, dy, dz);
            max_error_from_double = std::max(max_error_from_double, error);

            dx = geod_x_midpt - sphere_mid_x_f;
            dy = geod_y_midpt - sphere_mid_y_f;
            dz = geod_z_midpt - sphere_mid_z_f;

            error = length(dx, dy, dz);
            max_error_from_float = std::max(max_error_from_float, error);
        }

        printf("%4d  %8d %9.04f %9.04f\n",
            (int)line_divisor, (int)seg_len, max_error_from_double, max_error_from_float);
    }
}

void checkDoubleVsFloat(double lat1, double lon1, double lat2, double lon2)
{
    const double a = Constants::WGS84_a();
    const double f = Constants::WGS84_f();
    const Geodesic geods(a, f);
    const Geocentric earth(a, f);

    const double EARTH_CIRCUMFERENCE = 2.0 * GeographicLib::Math::pi() * a;

    double total_dist, start_dir, end_dir;
    geods.Inverse(lat1, lon1, lat2, lon2, total_dist, start_dir, end_dir);

    const double LINE_DIVISOR = 512;

    double max_seg_len = (EARTH_CIRCUMFERENCE / 2.0) / LINE_DIVISOR;
    int seg_count = (int)ceil(total_dist / max_seg_len);
    double seg_len = total_dist / (double)seg_count;

    double max_error = 0;
    for (int pt = 1; pt < seg_count; pt++)
    {
        GeoCoord ctrl_pt_d;
        geods.Direct(lat1, lon1, start_dir, seg_len * pt, ctrl_pt_d.lat, ctrl_pt_d.lon);

        GeoCoord ctrl_pt_f;
        geods.Direct(
            (float)lat1, (float)lon1, (float)start_dir, (float)seg_len * (float)pt,
            ctrl_pt_f.lat, ctrl_pt_f.lon);
        ctrl_pt_f.lat = (float)ctrl_pt_f.lat;
        ctrl_pt_f.lon = (float)ctrl_pt_f.lon;

        double err_dist, err_start_dir, err_end_dir;
        geods.Inverse(
            ctrl_pt_d.lat, ctrl_pt_d.lon, ctrl_pt_f.lat, ctrl_pt_f.lon,
            err_dist, err_start_dir, err_end_dir);
        max_error = std::max(max_error, err_dist);
    }

    printf("max float error: %f\n", max_error);
}

void check_coords(double lat1, double lon1, double lat2, double lon2)
{
    printf("\n");
    printf("(%f, %f) -> (%f, %f)\n", lat1, lon1, lat2, lon2);
    checkGeoVsSphere(lat1, lon1, lat2, lon2);
    checkDoubleVsFloat(lat1, lon1, lat2, lon2);
}

int _tmain(int argc, _TCHAR* argv[])
{
    double lat1 = 0;
    double lon1 = 0;
    double lat2 = 0;
    double lon2 = 179.9;
    check_coords(lat1, lon1, lat2, lon2);

    srand(0);

    for (int i = 0; i < 5; i++)
    {
        lat1 = ((double)rand() / (double)RAND_MAX) * 180.0 - 90.0;
        lon1 = ((double)rand() / (double)RAND_MAX) * 360.0 - 180.0;
        lat2 = ((double)rand() / (double)RAND_MAX) * 180.0 - 90.0;
        lon2 = ((double)rand() / (double)RAND_MAX) * 360.0 - 180.0;
        check_coords(lat1, lon1, lat2, lon2);
    }

    return 0;
}

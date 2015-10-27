#include "stdafx.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>

#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/Rhumb.hpp>

using namespace GeographicLib;

#define WGS84_R2_METERS 6371007.1809   // radius of sphere of equal area

template<typename T>
T metersToRad(T meters) {
    return (T)meters *  ((T)1.0 / (T)WGS84_R2_METERS);
}

template<typename T>
struct GeodeticCoord {
    T lat, lon, alt;
};

template<typename T>
struct GeocentricCoord {
    T x, y, z;

    GeocentricCoord<T> operator+(const GeocentricCoord<T>& rhs) const {
        return GeocentricCoord<T>({
            x + rhs.x,
            y + rhs.y,
            z + rhs.z
        });
    }

    GeocentricCoord<T> operator-(const GeocentricCoord<T>& rhs) const {
        return GeocentricCoord<T>({
            x - rhs.x,
            y - rhs.y,
            z - rhs.z
        });
    }

    GeocentricCoord<T> operator*(T scale) const {
        return GeocentricCoord<T>({
            x * scale,
            y * scale,
            z * scale
        });
    }

    GeocentricCoord<T> operator/(T scale) const {
        return GeocentricCoord<T>({
            x / scale,
            y / scale,
            z / scale
        });
    }
};

template<typename T>
struct Mat3x3 {
    T m11, m21, m31;
    T m12, m22, m32;
    T m13, m23, m33;

    GeocentricCoord<T> operator*(const GeocentricCoord<T>& coord) const {
        return GeocentricCoord<T>({
            m11 * coord.x + m21 * coord.y + m31 * coord.z,
            m12 * coord.x + m22 * coord.y + m32 * coord.z,
            m13 * coord.x + m23 * coord.y + m33 * coord.z
        });
    }
};

const double a = Constants::WGS84_a();
const double f = Constants::WGS84_f();
const double b = (1.0 - f) * a;
const Geodesic geods(a, f);
const Rhumb rhumb(a, f);
const Geocentric earth(a, f);
const double EARTH_CIRCUMFERENCE = 2.0 * GeographicLib::Math::pi() * a;
const double HALF_EARTH_CIRCUMFERENCE = EARTH_CIRCUMFERENCE * 0.5;

template<typename T>
T length(T x, T y, T z) {
    return sqrt(x * x + y * y + z * z);
}

template<typename T>
T length(GeocentricCoord<T> v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

template<typename T>
void normalize(T x, T y, T z, T& nx, T& ny, T& nz) {
    T len = length(x, y, z);
    nx = x / len;
    ny = y / len;
    nz = z / len;
}

template<typename T>
GeocentricCoord<T> normalize(GeocentricCoord<T> v) {
    return v / length(v);
}

template<typename T>
GeocentricCoord<T> cross(GeocentricCoord<T> v1, GeocentricCoord<T> v2) {
    GeocentricCoord<T> cross = {
        v1.y * v2.z - v1.z * v2.y,
        v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x
    };
    return cross;
}

template<typename T>
void calcRhumbLineDirectFvMethod(
    T lat1, T lon1,
    T distance, T heading,
    T& lat2, T& lon2)
{
    T rad_lat1, rad_lon1;
    T rad_lat2, rad_lon2;
    T rad_distance;
    T rad_heading;
    T dY, dX, q;

    // Fudge the starting point at the poles.
    lat1 = std::min(lat1, (T)89.999999);
    lat1 = std::max(lat1, (T)-89.999999);

    // convert point 1 to radians
    rad_lat1 = lat1 * (T)Math::degree();
    rad_lon1 = lon1 * (T)Math::degree();

    // input is in meters and degrees
    rad_distance = metersToRad(distance);
    rad_heading = heading * (T)Math::degree();

    // compute ending latitude
    rad_lat2 = rad_lat1 + rad_distance * cos(rad_heading);

    // Compute dY in a mercator projection with a scale factor of 
    // 1 / WGS84_R2_METERS.
    dY = log(tan(rad_lat2 / 2 + (T)Math::pi() / 4) / tan(rad_lat1 / 2 + (T)Math::pi() / 4));

    // Compute the inverse of the mercator projection stretching factor.
    if (fabs(rad_lat2 - rad_lat1) < sqrt(1.0e-15))
        q = cos(rad_lat1);
    else
        q = (rad_lat2 - rad_lat1) / dY;

    // Compute dX in a mercator projection with a scale factor of 
    // 1 / WGS84_R2_METERS.
    dX = rad_distance * sin(rad_heading) / q;

    // Compute ending longitude.  Must handle IDL crossing.
    rad_lon2 = rad_lon1 + dX;
    if (rad_lon2 > (T)Math::pi())
        rad_lon2 -= (T)Math::pi() * (T)2;
    else if (rad_lon2 < (T)-Math::pi())
        rad_lon2 += (T)Math::pi() * (T)2;

    lat2 = rad_lat2 / (T)Math::degree();
    lon2 = rad_lon2 / (T)Math::degree();
}

template<typename T>
GeocentricCoord<T> geoToXYZ(T lat, T lon, T alt)
{
    T n, a, b, e;
    T sinlat, coslat, sinlon, coslon;

    // Convert to radians
    T rlat = lat * (T)0.0174532925199433;
    T rlon = lon * (T)0.0174532925199433;

    sinlat = sin(rlat);
    coslat = cos(rlat);
    sinlon = sin(rlon);
    coslon = cos(rlon);

    a = 6378137.0;  // WGS84 meters at equator

    b = (T)6356752.3142;
    e = (T)0.99330561999;
    n = (a * a) / sqrt((a * a * coslat * coslat) + (b * b * sinlat * sinlat));

    GeocentricCoord<T> result;
    result.x = (n + alt) * coslat * coslon;
    result.y = (n + alt) * coslat * sinlon;
    result.z = ((n * e) + alt) * sinlat;

    return result;
}

template<typename T>
void getMidpointXyz(GeodeticCoord<double> pt1, GeodeticCoord<double> pt2,
    T& sphere_mid_x, T& sphere_mid_y, T& sphere_mid_z) {

    double x1_d, y1_d, z1_d, x2_d, y2_d, z2_d;
    earth.Forward(pt1.lat, pt1.lon, pt1.alt, x1_d, y1_d, z1_d);
    earth.Forward(pt2.lat, pt2.lon, pt2.alt, x2_d, y2_d, z2_d);

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

double getMidpointDelta(GeodeticCoord<double> pt1, GeodeticCoord<double> pt2)
{
    GeocentricCoord<double> p1, p2;
    earth.Forward(pt1.lat, pt1.lon, pt1.alt, p1.x, p1.y, p1.z);
    earth.Forward(pt2.lat, pt2.lon, pt2.alt, p2.x, p2.y, p2.z);

    double dist, start_dir, end_dir;
    geods.Inverse(pt1.lat, pt1.lon, pt2.lat, pt2.lon, dist, start_dir, end_dir);

    double lat2, lon2;
    geods.Direct(pt1.lat, pt1.lon, start_dir, dist * 0.5, lat2, lon2);

    GeocentricCoord<double> mid_geod;
    earth.Forward(lat2, lon2, (pt1.alt + pt2.alt) * 0.5, mid_geod.x, mid_geod.y, mid_geod.z);

    GeocentricCoord<double> mid_straight = (p1 + p2) * 0.5;

    return length(mid_geod - mid_straight);
}

enum LineType { GC, RHUMB };

void getError(
    double lat1, double lon1, double lat2, double lon2, LineType line_type,
    double total_dist, double start_dir, double line_divisor,
    double& seg_len, double& max_error) {

    const double ALTITUDE = 0;

    double max_seg_len = (EARTH_CIRCUMFERENCE / 2.0) / line_divisor;
    int seg_count = (int)ceil(total_dist / max_seg_len);
    seg_len = total_dist / (double)seg_count;

    std::vector<GeodeticCoord<double>> coords;
    coords.push_back({ lat1, lon1, ALTITUDE });

    for (int pt = 1; pt < seg_count; pt++) {
        GeodeticCoord<double> ctrl_pt;
        if (line_type == GC) {
            geods.Direct(lat1, lon1, start_dir, seg_len * pt, ctrl_pt.lat, ctrl_pt.lon);
        }
        else {
            rhumb.Direct(lat1, lon1, start_dir, seg_len * pt, ctrl_pt.lat, ctrl_pt.lon);
        }
        ctrl_pt.alt = ALTITUDE;
        coords.push_back(ctrl_pt);
    }

    coords.push_back({ lat2, lon2, ALTITUDE });

    max_error = 0;
    for (unsigned int ctrl_pt = 0; ctrl_pt < coords.size() - 1; ctrl_pt++) {
        double geod_dist, geod_start_dir, geod_end_dir;
        auto& pt1 = coords[ctrl_pt];
        auto& pt2 = coords[ctrl_pt + 1];
        if (line_type == GC) {
            geods.Inverse(
                pt1.lat, pt1.lon, pt2.lat, pt2.lon,
                geod_dist, geod_start_dir, geod_end_dir);
        }
        else {
            rhumb.Inverse(
                pt1.lat, pt1.lon, pt2.lat, pt2.lon,
                geod_dist, geod_start_dir, geod_end_dir);
        }

        GeodeticCoord<double> geod_latlon_midpt;
        if (line_type == GC) {
            geods.Direct(
                pt1.lat, pt1.lon, geod_start_dir, geod_dist / 2.0,
                geod_latlon_midpt.lat, geod_latlon_midpt.lon);
        }
        else {
            rhumb.Direct(
                pt1.lat, pt1.lon, geod_start_dir, geod_dist / 2.0,
                geod_latlon_midpt.lat, geod_latlon_midpt.lon);
        }

        geod_latlon_midpt.alt = (pt1.alt + pt2.alt) / 2.0;

        double geod_x_midpt, geod_y_midpt, geod_z_midpt;
        earth.Forward(
            geod_latlon_midpt.lat, geod_latlon_midpt.lon, geod_latlon_midpt.alt,
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
        max_error = std::max(max_error, error);

        //dx = geod_x_midpt - sphere_mid_x_f;
        //dy = geod_y_midpt - sphere_mid_y_f;
        //dz = geod_z_midpt - sphere_mid_z_f;

        //error = length(dx, dy, dz);
        //max_error_from_float = std::max(max_error_from_float, error);
    }
}

void checkGeoVsSphere(double lat1, double lon1, double lat2, double lon2) {
    double total_dist_gc, start_dir_gc, end_dir_gc;
    geods.Inverse(lat1, lon1, lat2, lon2, total_dist_gc, start_dir_gc, end_dir_gc);

    double total_dist_rhumb, start_dir_rhumb, end_dir_rhumb;
    rhumb.Inverse(lat1, lon1, lat2, lon2, total_dist_rhumb, start_dir_rhumb, end_dir_rhumb);

    printf("divs seg dist(g) seg dist(r)   err (g)      err (r)\n");
    printf("---- ----------- ----------- --------- ------------\n");

    for (double line_divisor = 2;
        line_divisor <= 4096;
        line_divisor *= 2) {
        double seg_len_gc, seg_len_rhumb, max_error_from_gc, max_error_from_rhumb;

        getError(
            lat1, lon1, lat2, lon2, GC,
            total_dist_gc, start_dir_gc, line_divisor,
            seg_len_gc, max_error_from_gc);

        getError(
            lat1, lon1, lat2, lon2, RHUMB,
            total_dist_rhumb, start_dir_rhumb, line_divisor,
            seg_len_rhumb, max_error_from_rhumb);

        printf("%4d %11d %11d %9.04f %12.04f\n",
            (int)line_divisor, (int)seg_len_gc, (int)seg_len_rhumb,
            max_error_from_gc, max_error_from_rhumb);
    }
}

void checkDoubleVsFloat(double lat1, double lon1, double lat2, double lon2) {
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
    for (int pt = 1; pt < seg_count; pt++) {
        GeodeticCoord<double> ctrl_pt_d;
        geods.Direct(lat1, lon1, start_dir, seg_len * pt, ctrl_pt_d.lat, ctrl_pt_d.lon);

        GeodeticCoord<double> ctrl_pt_f;
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

void check_coords(double lat1, double lon1, double lat2, double lon2) {
    printf("\n");
    printf("(%f, %f) -> (%f, %f)\n", lat1, lon1, lat2, lon2);
    checkGeoVsSphere(lat1, lon1, lat2, lon2);
    checkDoubleVsFloat(lat1, lon1, lat2, lon2);
}

template<typename T>
void generateArcTemplate(
    int line_divisor, int seg_divisor,
    std::vector<GeocentricCoord<T>>& arc_template) {

    GeocentricCoord<T> coord = { 0, 0, 0 };
    arc_template.clear();

    const T ARC_LEN = (T)2.0 * (T)M_PI * (T)0.5 * ((T)1.0 / ((T)line_divisor * (T)seg_divisor));
    for (int i = 0; i <= seg_divisor; i++) {
        T angle = (T)i * ARC_LEN;
        coord.x = cos(angle) - (T)1; // subtract 1 to shift template so p1 is at origin
        coord.y = sin(angle);
        arc_template.push_back(coord);
    }
}

template<typename T>
Mat3x3<T> getTempalteToSegMatrix(
    GeocentricCoord<T> p1, GeocentricCoord<T> p2, GeocentricCoord<T> p1_nrm)
{
    GeocentricCoord<T> seg_x_axis = p1_nrm;
    GeocentricCoord<T> seg_y_axis = p2 - p1;
    GeocentricCoord<T> seg_z_axis = cross(seg_x_axis, seg_y_axis);
    seg_z_axis = normalize(seg_z_axis);
    seg_y_axis = cross(seg_z_axis, seg_x_axis);

    Mat3x3<T> template_to_seg = {
        seg_x_axis.x, seg_y_axis.x, seg_z_axis.x,
        seg_x_axis.y, seg_y_axis.y, seg_z_axis.y,
        seg_x_axis.z, seg_y_axis.z, seg_z_axis.z
    };

    return template_to_seg;
}

template<typename T>
GeocentricCoord<T> getCoordFromTemplate(
    GeocentricCoord<T> template_coord,
    T scale_factor, Mat3x3<T> template_to_seg_transform,
    GeocentricCoord<T> p1) {

    template_coord = template_coord * scale_factor;
    template_coord = template_to_seg_transform * template_coord;
    template_coord = template_coord + p1;

    return template_coord;
}


void testArcTemplate(
    GeodeticCoord<double> coord1, GeodeticCoord<double> coord2,
    const std::vector<GeocentricCoord<double>>& arc_template_d,
    const std::vector<GeocentricCoord<float>>& arc_template_f,
    double& max_error_via_f_gc,
    //double& max_error_via_f_ut_gc,
    double& max_error_via_d_gc,
    double& max_midpoint_error_d_gc,
    double& max_error_via_f_rhumb,
    //double& max_error_via_f_ut_rhumb,
    double& max_error_via_d_rhumb) {

    //printf("(%7.02f, %7.02f, %d)->(%7.02f, %7.02f, %d)\n",
    //    coord1.lat, coord1.lon, (int)coord1.alt,
    //    coord2.lat, coord2.lon, (int)coord2.alt);

    double base_dist_gc, start_dir_gc, end_dir_gc;
    geods.Inverse(
        coord1.lat, coord1.lon, coord2.lat, coord2.lon,
        base_dist_gc, start_dir_gc, end_dir_gc);

   double base_dist_rhumb, start_dir_rhumb, end_dir_rhumb;
    rhumb.Inverse(
        coord1.lat, coord1.lon, coord2.lat, coord2.lon,
        base_dist_rhumb, start_dir_rhumb, end_dir_rhumb);

    const int SUBSEG_COUNT = arc_template_f.size() - 1;

    std::vector<GeocentricCoord<double>> baseline_wpts_gc;
    std::vector<GeocentricCoord<double>> baseline_wpts_rhumb;
    for (int i = 0; i <= SUBSEG_COUNT; i++) {
        double dist;
        double lat, lon, alt;
        GeocentricCoord<double> wpt;

        double frac = (double)i / (double)SUBSEG_COUNT;

        // gc baseline wpt
        dist = base_dist_gc * frac;

        geods.Direct(coord1.lat, coord1.lon, start_dir_gc, dist, lat, lon);
        alt = coord1.alt + (coord2.alt - coord1.alt) * frac;

        earth.Forward(lat, lon, alt, wpt.x, wpt.y, wpt.z);
        baseline_wpts_gc.push_back(wpt);

        // rhumb baseline wpt
        dist = base_dist_rhumb * frac;

        geods.Direct(coord1.lat, coord1.lon, start_dir_rhumb, dist, lat, lon);
        alt = coord1.alt + (coord2.alt - coord1.alt) * frac;

        earth.Forward(lat, lon, alt, wpt.x, wpt.y, wpt.z);
        baseline_wpts_rhumb.push_back(wpt);
    }

    GeocentricCoord<double> p1_d, p2_d;
    earth.Forward(coord1.lat, coord1.lon, coord1.alt, p1_d.x, p1_d.y, p1_d.z);
    earth.Forward(coord2.lat, coord2.lon, coord2.alt, p2_d.x, p2_d.y, p2_d.z);

    const double EYE_DIST = 160000.0;
    GeocentricCoord<double> eye_pos = { p1_d.x + EYE_DIST, p1_d.y, p1_d.z };

    GeocentricCoord<float> p1_f = {
        (float)(p1_d.x - eye_pos.x),
        (float)(p1_d.y - eye_pos.y),
        (float)(p1_d.z - eye_pos.z)
    };
    GeocentricCoord<float> p2_f = {
        (float)(p2_d.x - eye_pos.x),
        (float)(p2_d.y - eye_pos.y),
        (float)(p2_d.z - eye_pos.z)
    };

    GeocentricCoord<float> p1_f_ut = {
        (float)p1_d.x,
        (float)p1_d.y,
        (float)p1_d.z
    };
    GeocentricCoord<float> p2_f_ut = {
        (float)p2_d.x,
        (float)p2_d.y,
        (float)p2_d.z
    };

    GeocentricCoord<double> p1_nrm_d;
    earth.Forward(coord1.lat, coord1.lon, coord1.alt + 10, p1_nrm_d.x, p1_nrm_d.y, p1_nrm_d.z);
    p1_nrm_d.x -= p1_d.x;
    p1_nrm_d.y -= p1_d.y;
    p1_nrm_d.z -= p1_d.z;
    p1_nrm_d = normalize(p1_nrm_d);

    GeocentricCoord<float> p1_nrm_f = {
        (float)p1_nrm_d.x,
        (float)p1_nrm_d.y,
        (float)p1_nrm_d.z
    };

    p1_d = p1_d - eye_pos;
    p2_d = p2_d - eye_pos;

    Mat3x3<float> template_to_seg_f = getTempalteToSegMatrix(p1_f, p2_f, p1_nrm_f);
    Mat3x3<double> template_to_seg_d = getTempalteToSegMatrix(p1_d, p2_d, p1_nrm_d);

    max_error_via_f_gc = 0;
    //max_error_via_f_ut_gc = 0;
    max_error_via_d_gc = 0;
    max_midpoint_error_d_gc = 0;

    float template_scale_factor_f = length(p2_f - p1_f) / length(arc_template_f.back());
    //float template_scale_factor_f_ut = length(p2_f_ut - p1_f_ut) / length(arc_template_f.back());
    double template_scale_factor_d = length(p2_d - p1_d) / length(arc_template_d.back());

    GeocentricCoord<float> template_coord_back_f = arc_template_f.back();
    template_coord_back_f = getCoordFromTemplate(
        template_coord_back_f, template_scale_factor_f, template_to_seg_f, p1_f);

    //GeocentricCoord<float> template_coord_back_f_ut = arc_template_f.back();
    //template_coord_back_f_ut = getCoordFromTemplate(
    //    template_coord_back_f_ut, template_scale_factor_f_ut, template_to_seg_f, p1_f_ut);

    GeocentricCoord<double> template_coord_back_d = arc_template_d.back();
    template_coord_back_d = getCoordFromTemplate(
        template_coord_back_d, template_scale_factor_d, template_to_seg_d, p1_d);

    GeocentricCoord<float> offset_f = p2_f - template_coord_back_f;
    //GeocentricCoord<float> offset_f_ut = p2_f_ut - template_coord_back_f_ut;
    GeocentricCoord<double> offset_d = p2_d - template_coord_back_d;

    GeocentricCoord<double> last_template_coord_d;

    for (unsigned int i = 0; i < baseline_wpts_gc.size(); i++) {
        float offset_frac_f = (float)i / ((float)baseline_wpts_gc.size() - 1.0f);
        double offset_frac_d = (double)i / ((double)baseline_wpts_gc.size() - 1.0);

        offset_frac_f *= offset_frac_f;
        offset_frac_d *= offset_frac_d;

        GeocentricCoord<float> template_coord_f = arc_template_f[i];
        template_coord_f = getCoordFromTemplate(
            template_coord_f, template_scale_factor_f, template_to_seg_f, p1_f);
        GeocentricCoord<float> correction_f = offset_f * offset_frac_f;
        template_coord_f = template_coord_f + correction_f;

        //GeocentricCoord<float> template_coord_f_ut = arc_template_f[i];
        //template_coord_f_ut = getCoordFromTemplate(
        //    template_coord_f_ut, template_scale_factor_f_ut, template_to_seg_f, p1_f_ut);
        //GeocentricCoord<float> correction_f_ut = offset_f_ut * offset_frac_f;
        //template_coord_f_ut = template_coord_f_ut + correction_f_ut;

        GeocentricCoord<double> template_coord_d = arc_template_d[i];
        template_coord_d = getCoordFromTemplate(
            template_coord_d, template_scale_factor_d, template_to_seg_d, p1_d);
        GeocentricCoord<double> correction_d = offset_d * offset_frac_d;
        template_coord_d = template_coord_d + correction_d;

        GeocentricCoord<double> template_coord_f_as_d = {
            (double)template_coord_f.x,
            (double)template_coord_f.y,
            (double)template_coord_f.z
        };

        //GeocentricCoord<double> template_coord_f_ut_as_d = {
        //    (double)template_coord_f_ut.x,
        //    (double)template_coord_f_ut.y,
        //    (double)template_coord_f_ut.z
        //};

        // get gc error
        GeocentricCoord<double> ref_coord_gc = baseline_wpts_gc[i];
        ref_coord_gc = ref_coord_gc - eye_pos;

        GeocentricCoord<double> ref_coord_ut_gc = baseline_wpts_gc[i];

        double error_via_f_gc = length(template_coord_f_as_d - ref_coord_gc);
        double error_via_d_gc = length(template_coord_d - ref_coord_gc);

        //double error_via_f_ut_gc = length(template_coord_f_ut_as_d - ref_coord_ut_gc);

        // get midpoint error
        if (i > 0) {
            GeocentricCoord<double> midpt = (last_template_coord_d + template_coord_d) * 0.5;

            double lat1, lon1, alt1;
            earth.Reverse(
                last_template_coord_d.x, last_template_coord_d.y, last_template_coord_d.z,
                lat1, lon1, alt1);

            double lat2, lon2, alt2;
            earth.Reverse(
                template_coord_d.x, template_coord_d.y, template_coord_d.z,
                lat2, lon2, alt2);

            double dist, start_dir, end_dir;
            geods.Inverse(lat1, lon1, lat2, lon2, dist, start_dir, end_dir);

            double mid_lat, mid_lon;
            geods.Direct(lat1, lon1, start_dir, dist * 0.5, mid_lat, mid_lon);

            GeocentricCoord<double> geod_midpt;
            earth.Forward(mid_lat, mid_lon, (alt1 + alt2) * 0.5,
                geod_midpt.x, geod_midpt.y, geod_midpt.z);

            double midpoint_error = length(geod_midpt - midpt);
            max_midpoint_error_d_gc = std::max(max_midpoint_error_d_gc, midpoint_error);
        }

        last_template_coord_d = template_coord_d;

        // get rhumb line error
        GeocentricCoord<double> ref_coord_rhumb = baseline_wpts_rhumb[i];
        ref_coord_rhumb = ref_coord_rhumb - eye_pos;

        //GeocentricCoord<double> ref_coord_ut_rhumb = baseline_wpts_rhumb[i];

        double error_via_f_rhumb = length(template_coord_f_as_d - ref_coord_rhumb);
        double error_via_d_rhumb = length(template_coord_d - ref_coord_rhumb);

        //double error_via_f_ut_rhumb = length(template_coord_f_ut_as_d - ref_coord_ut_rhumb);

        max_error_via_f_gc = std::max(max_error_via_f_gc, error_via_f_gc);
        //max_error_via_f_ut_gc = std::max(max_error_via_f_ut_gc, error_via_f_ut_gc);
        max_error_via_d_gc = std::max(max_error_via_d_gc, error_via_d_gc);

        max_error_via_f_rhumb = std::max(max_error_via_f_rhumb, error_via_f_rhumb);
        //max_error_via_f_ut_rhumb = std::max(max_error_via_f_ut_rhumb, error_via_f_ut_rhumb);
        max_error_via_d_rhumb = std::max(max_error_via_d_rhumb, error_via_d_rhumb);

        //printf("%f / %f\n", error_via_f, error_via_d);
        //printf("%f / %f / %f\n", error_via_f, error_via_f_ut, error_via_d);
        //printf("%f / %f / %f\n", error_via_f_rhumb, error_via_f_ut_rhumb, error_via_d_rhumb);
    }
    //printf("%f / %f\n", max_error_via_f, max_error_via_d);
}

void testArcTemplateMethod() {
    const int TEST_COUNT = 25;
    const int DIVISOR_PAIRS[][2] = {
        { 64, 64 },
        { 128, 32 },
        { 256, 16 },
        { 256, 64 },
        { 512, 8 }
    };

    const int PAIR_IDX = 2;
    const int STARTING_LINE_DIVISOR = DIVISOR_PAIRS[PAIR_IDX][0];
    const int SEG_DIVISOR = DIVISOR_PAIRS[PAIR_IDX][1];

    const double STARTING_SEG_LEN = HALF_EARTH_CIRCUMFERENCE / STARTING_LINE_DIVISOR;

    printf("max seg len: %f\n", STARTING_SEG_LEN);

    std::vector<GeocentricCoord<double>> arc_template_d;
    std::vector<GeocentricCoord<float>> arc_template_f;
    generateArcTemplate(STARTING_LINE_DIVISOR, SEG_DIVISOR, arc_template_d);
    generateArcTemplate(STARTING_LINE_DIVISOR, SEG_DIVISOR, arc_template_f);

    double max_error_via_f_gc = 0;
    //double max_error_via_f_ut_gc = 0;
    double max_error_via_d_gc = 0;
    double max_midpoint_error_d_gc = 0;
    double max_error_via_f_rhumb = 0;
    //double max_error_via_f_ut_rhumb = 0;
    double max_error_via_d_rhumb = 0;

    srand(0);
    for (int i = 0; i < TEST_COUNT; i++) {
        double dir = ((double)rand() / (double)RAND_MAX) * 360.0;
        double lat1 = ((double)rand() / (double)RAND_MAX) * 180.0 - 90.0;
        double lon1 = ((double)rand() / (double)RAND_MAX) * 360.0 - 180.0;

        double lat2, lon2;
        geods.Direct(lat1, lon1, dir, STARTING_SEG_LEN, lat2, lon2);

        double alt1 = ((double)rand() / (double)RAND_MAX) * 100000.0;
        double alt2 = alt1;

        double max_error_via_f_inner_gc = 0;
        //double max_error_via_f_ut_inner_gc = 0;
        double max_error_via_d_inner_gc = 0;
        double max_midpoint_error_d_inner_gc = 0;
        double max_error_via_f_inner_rhumb = 0;
        //double max_error_via_f_ut_inner_rhumb = 0;
        double max_error_via_d_inner_rhumb = 0;
        testArcTemplate(
            { lat1, lon1, alt1 }, { lat2, lon2, alt2 },
            arc_template_d, arc_template_f,
            max_error_via_f_inner_gc,
            //max_error_via_f_ut_inner_gc,
            max_error_via_d_inner_gc,
            max_midpoint_error_d_inner_gc,
            max_error_via_f_inner_rhumb,
            //max_error_via_f_ut_inner_rhumb,
            max_error_via_d_inner_rhumb);

        max_error_via_f_gc = std::max(max_error_via_f_gc, max_error_via_f_inner_gc);
        //max_error_via_f_ut_gc = std::max(max_error_via_f_ut_gc, max_error_via_f_ut_inner_gc);
        max_error_via_d_gc = std::max(max_error_via_d_gc, max_error_via_d_inner_gc);
        max_midpoint_error_d_gc = std::max(max_midpoint_error_d_gc, max_midpoint_error_d_inner_gc);
        max_error_via_f_rhumb = std::max(max_error_via_f_rhumb, max_error_via_f_inner_rhumb);
        //max_error_via_f_ut_rhumb = std::max(max_error_via_f_ut_rhumb, max_error_via_f_ut_inner_rhumb);
        max_error_via_d_rhumb = std::max(max_error_via_d_rhumb, max_error_via_d_inner_rhumb);
    }

    //printf("\n%.3f / %.3f / %.3f / %.3f / %.3f / %.3f / %.3f\n",
    //    max_error_via_f_gc, max_error_via_f_ut_gc, max_error_via_d_gc, max_midpoint_error_d_gc,
    //    max_error_via_f_rhumb, max_error_via_f_ut_rhumb, max_error_via_d_rhumb);
    printf("\n%.3f / %.3f / %.3f / %.3f / %.3f\n",
        max_error_via_f_gc, max_error_via_d_gc, max_midpoint_error_d_gc,
        max_error_via_f_rhumb, max_error_via_d_rhumb);
}

void testGeodMidpointToStraightLineMidpointError() {
    const int TEST_COUNT = 100000;
    const int STARTING_LINE_DIVISOR = 256;
    const double STARTING_SEG_LEN = HALF_EARTH_CIRCUMFERENCE / STARTING_LINE_DIVISOR;
    const int NUM_LENGTHS_TO_TEST = 8;

    double max_errors[NUM_LENGTHS_TO_TEST];
    memset(max_errors, 0, sizeof(max_errors));

    srand(0);
    for (int i = 0; i < TEST_COUNT; i++) {
        double dir = ((double)rand() / (double)RAND_MAX) * 360.0;
        double lat_mid = ((double)rand() / (double)RAND_MAX) * 180.0 - 90.0;
        double lon_mid = ((double)rand() / (double)RAND_MAX) * 360.0 - 180.0;
        //printf("----\n");
        for (int j = 0; j < NUM_LENGTHS_TO_TEST; j++) {
            double half_dist = STARTING_SEG_LEN / pow(2.0, (double)(j + 1));
            GeodeticCoord<double> p1, p2;
            geods.Direct(lat_mid, lon_mid, dir, -half_dist, p1.lat, p1.lon);
            geods.Direct(lat_mid, lon_mid, dir, half_dist, p2.lat, p2.lon);
            p1.alt = 100000;
            p2.alt = 100000;
            double delta = getMidpointDelta(p1, p2);
            max_errors[j] = std::max(max_errors[j], delta);
            //printf("delta for %d segs: %f\n", STARTING_LINE_DIVISOR * (int)pow(2, j), delta);
        }
    }

    printf("max deltas for for %d tests:\n", TEST_COUNT);
    for (int i = 0; i < NUM_LENGTHS_TO_TEST; i++) {
        double dist = STARTING_SEG_LEN / pow(2.0, (double)i);
        printf("%5d segs of length %5d m: %10f m\n",
            STARTING_LINE_DIVISOR * (int)pow(2, i), (int)round(dist), max_errors[i]);
    }
}

void testRhumbLineAsFloatError()
{
    const int TEST_COUNT = 100000;
    const int STARTING_LINE_DIVISOR = 64;
    const int SUBSEG_COUNT = 64;
    const double STARTING_SEG_LEN = HALF_EARTH_CIRCUMFERENCE / STARTING_LINE_DIVISOR;

    double max_err_fv_method_from_f = 0;
    double max_err_fv_method_from_d = 0;

    srand(0);
    for (int i = 0; i < TEST_COUNT; i++) {

        if ((i + 1) % 100 == 0) {
            printf("\r%d tests done...", i + 1);
        }

        GeodeticCoord<double> p1, p2;

        double dir = ((double)rand() / (double)RAND_MAX) * 360.0;
        p1.lat = ((double)rand() / (double)RAND_MAX) * 180.0 - 90.0;
        p1.lon = ((double)rand() / (double)RAND_MAX) * 360.0 - 180.0;
        p1.alt = ((double)rand() / (double)RAND_MAX) * 100000.0;

        rhumb.Direct(p1.lat, p1.lon, dir, STARTING_SEG_LEN, p2.lat, p2.lon);
        p2.alt = p1.alt;

        double max_err_fv_method_from_f_inner = 0;
        double max_err_fv_method_from_d_inner = 0;

        float dlat = (float)p2.lat - (float)p1.lat;
        float dlon = (float)p2.lon - (float)p1.lon;
        if (dlon > 180.0f) {
            dlon -= 360.0f;
        }
        if (dlon < -180.0f) {
            dlon += 360.0f;
        }
        float dalt = (float)p2.alt - (float)p1.alt;
        for (int j = 0; j <= SUBSEG_COUNT; j++) {
            GeodeticCoord<float> interp_f;
            float frac_f = (float)j / (float)SUBSEG_COUNT;
            float interp_dist_f = frac_f * (float)STARTING_SEG_LEN;

            calcRhumbLineDirectFvMethod(
                (float)p1.lat, (float)p1.lon, interp_dist_f, (float)dir,
                interp_f.lat, interp_f.lon);
            interp_f.alt = (float)p1.alt + frac_f * dalt;

            double frac_d = (double)j / (double)SUBSEG_COUNT;
            double interp_dist_d = frac_d * STARTING_SEG_LEN;

            GeodeticCoord<double> interp_d;
            calcRhumbLineDirectFvMethod(
                p1.lat, p1.lon, interp_dist_d, dir,
                interp_d.lat, interp_d.lon);
            interp_d.alt = p1.alt + frac_d * dalt;

            GeodeticCoord<double> interp_d_ref;
            rhumb.Direct(p1.lat, p1.lon, dir, interp_dist_d, interp_d_ref.lat, interp_d_ref.lon);
            interp_d_ref.alt = p1.alt + frac_d * dalt;

            GeocentricCoord<float> xyz_f = geoToXYZ(interp_f.lat, interp_f.lon, interp_f.alt);
            GeocentricCoord<double> xyz_d = geoToXYZ(interp_d.lat, interp_d.lon, interp_d.alt);

            GeocentricCoord<double> xyz_d_ref;
            earth.Forward(
                interp_d_ref.lat, interp_d_ref.lon, interp_d_ref.alt,
                xyz_d_ref.x, xyz_d_ref.y, xyz_d_ref.z);

            GeocentricCoord<double> xyz_f_as_d = { xyz_f.x, xyz_f.y, xyz_f.z };

            double err_fv_method_from_f = length(xyz_f_as_d - xyz_d_ref);
            double err_fv_method_from_d = length(xyz_d - xyz_d_ref);
            max_err_fv_method_from_f_inner =
                std::max(max_err_fv_method_from_f_inner, err_fv_method_from_f);
            max_err_fv_method_from_d_inner =
                std::max(max_err_fv_method_from_d_inner, err_fv_method_from_d);
        }

        max_err_fv_method_from_f =
            std::max(max_err_fv_method_from_f_inner, max_err_fv_method_from_f_inner);
        max_err_fv_method_from_d =
            std::max(max_err_fv_method_from_d_inner, max_err_fv_method_from_d_inner);

        //printf("err for (%7.3f, %8.3f)->(%7.3f, %8.3f): %8.3f / %8.3f\n",
        //    p1.lat, p1.lon, p2.lat, p2.lon,
        //    max_err_fv_method_from_f_inner, max_err_fv_method_from_d_inner);
    }

    printf("\nmax err from FV method over %d tests (float/double): %8.3f m / %8.3f m\n",
        TEST_COUNT, max_err_fv_method_from_f, max_err_fv_method_from_d);
}

void testGeoToXYZFloatVsDouble() {
    const int TEST_COUNT = 100000;

    double max_err = 0;

    srand(0);
    for (int i = 0; i < TEST_COUNT; i++) {

        GeodeticCoord<double> p1_d;
        p1_d.lat = ((double)rand() / (double)RAND_MAX) * 180.0 - 90.0;
        p1_d.lon = ((double)rand() / (double)RAND_MAX) * 360.0 - 180.0;
        p1_d.alt = ((double)rand() / (double)RAND_MAX) * 100000.0;

        GeodeticCoord<float> p1_f = {
            (float)p1_d.lat, (float)p1_d.lon, (float)p1_d.alt
        };

        GeocentricCoord<float> xyz_f = geoToXYZ(p1_f.lat, p1_f.lon, p1_f.alt);
        GeocentricCoord<double> xyz_d = geoToXYZ(p1_d.lat, p1_d.lon, p1_d.alt);

        GeocentricCoord<double> xyz_f_as_d = {
            xyz_f.x, xyz_f.y, xyz_f.z
        };

        double err = length(xyz_d - xyz_f_as_d);
        max_err = std::max(max_err, err);
        //printf("err: %f\n", err);
    }

    printf("max err over %d tests: %f m\n", TEST_COUNT, max_err);
}

int _tmain(int argc, _TCHAR* argv[]) {
    double lat1 = 0;
    double lon1 = 0;
    double lat2 = 0;
    double lon2 = 179.9;
    //check_coords(lat1, lon1, lat2, lon2);

    srand(0);

    //for (int i = 0; i < 5; i++)
    //{
    //    lat1 = ((double)rand() / (double)RAND_MAX) * 180.0 - 90.0;
    //    lon1 = ((double)rand() / (double)RAND_MAX) * 360.0 - 180.0;
    //    lat2 = ((double)rand() / (double)RAND_MAX) * 180.0 - 90.0;
    //    lon2 = ((double)rand() / (double)RAND_MAX) * 360.0 - 180.0;
    //    check_coords(lat1, lon1, lat2, lon2);
    //}





    //testArcTemplateMethod();

    //testGeodMidpointToStraightLineMidpointError();

   //testRhumbLineAsFloatError();

    testGeoToXYZFloatVsDouble();

    return 0;
}

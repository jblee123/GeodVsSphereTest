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

const int LINE_DIVISOR = 256;
const int SEG_DIVISOR = 16;

const double MAX_SEG_LEN = EARTH_CIRCUMFERENCE * 0.5 / LINE_DIVISOR;

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
void generateArcTemplate(std::vector<GeocentricCoord<T>>& arc_template) {
	GeocentricCoord<T> coord = { 0, 0, 0 };
	arc_template.clear();

	const T ARC_LEN = (T)2.0 * (T)M_PI * (T)0.5 * ((T)1.0 / ((T)LINE_DIVISOR * (T)SEG_DIVISOR));
	for (int i = 0; i <= 16; i++) {
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

void testSegment(GeodeticCoord<double> coord1, GeodeticCoord<double> coord2,
	const std::vector<GeocentricCoord<double>>& arc_template_d,
	const std::vector<GeocentricCoord<float>>& arc_template_f) {
	printf("(%f, %f, %f)->(%f, %f, %f)\n",
		coord1.lat, coord1.lon, coord1.alt,
		coord2.lat, coord2.lon, coord2.alt);

	double base_dist, start_dir, end_dir;
	geods.Inverse(coord1.lat, coord1.lon, coord2.lat, coord2.lon, base_dist, start_dir, end_dir);

	std::vector<GeocentricCoord<double>> baseline_wpts;
	for (int i = 1; i <= SEG_DIVISOR - 1; i++) {
		double frac = (double)i / (double)SEG_DIVISOR;
		double dist = base_dist * frac;

		double lat, lon, alt;
		geods.Direct(coord1.lat, coord1.lon, start_dir, dist, lat, lon);
		alt = coord1.alt + (coord2.alt - coord1.alt) * frac;

		GeocentricCoord<double> wpt;
		earth.Forward(lat, lon, alt, wpt.x, wpt.y, wpt.z);
		baseline_wpts.push_back(wpt);
	}

	//const float EXPAND_FACTOR = (float)a / (float)b;
	//const float SHRINK_FACTOR = (float)b / (float)a;

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

	float template_scale_factor_f = length(p2_f - p1_f) / length(arc_template_f.back());
	double template_scale_factor_d = length(p2_d - p1_d) / length(arc_template_d.back());
	for (int i = 0; i < baseline_wpts.size(); i++) {
		GeocentricCoord<float> template_coord_f = arc_template_f[i + 1];
		template_coord_f = template_coord_f * template_scale_factor_f;
		template_coord_f = template_to_seg_f * template_coord_f;
		template_coord_f = template_coord_f + p1_f;

		GeocentricCoord<double> template_coord_d = arc_template_d[i + 1];
		template_coord_d = template_coord_d * template_scale_factor_d;
		template_coord_d = template_to_seg_d * template_coord_d;
		template_coord_d = template_coord_d + p1_d;

		GeocentricCoord<double> template_coord_f_as_d = {
			(double)template_coord_f.x,
			(double)template_coord_f.y,
			(double)template_coord_f.z
		};

		GeocentricCoord<double> ref_coord = baseline_wpts[i];
		ref_coord = ref_coord - eye_pos;

		double error_via_f = length(template_coord_f_as_d - ref_coord);
		double error_via_d = length(template_coord_d - ref_coord);
		printf("%f / %f\n", error_via_f, error_via_d);
	}
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

	std::vector<GeocentricCoord<double>> arc_template_d;
	std::vector<GeocentricCoord<float>> arc_template_f;
	generateArcTemplate(arc_template_d);
	generateArcTemplate(arc_template_f);

	for (int i = 0; i < 25; i++) {
		double dir = ((double)rand() / (double)RAND_MAX) * 360.0;
		lat1 = ((double)rand() / (double)RAND_MAX) * 180.0 - 90.0;
		lon1 = ((double)rand() / (double)RAND_MAX) * 360.0 - 180.0;
		geods.Direct(lat1, lon1, dir, MAX_SEG_LEN, lat2, lon2);

		double alt1 = ((double)rand() / (double)RAND_MAX) * 100000.0;
		double alt2 = alt1;

		//double alt1 = 0;
		//double alt2 = 0;
		//lat1 = 0;
		//lon1 = 0;
		//double dir = 45;
		//geods.Direct(lat1, lon1, dir, MAX_SEG_LEN, lat2, lon2);

		testSegment(
			{ lat1, lon1, alt1 }, { lat2, lon2, alt2 },
			arc_template_d, arc_template_f);
	}
	
	return 0;
}

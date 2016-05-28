#include "RhumbLines.h"

#include <stdbool.h>
#include <float.h>
#include <math.h>

#ifdef _WIN32
#define inline __inline
#endif

#define PI               3.14159265358979323846
#define DEG_TO_RAD(degrees)       (((double)(degrees)) * 1.7453292519943295e-2)
#define RAD_TO_DEG(radians)       (((double)(radians)) * 57.295779513082322)
#define WGS84_a_METERS  6378137.0      // Radius at the equator (semi-major)
#define WGS84_b_METERS  6356752.3142   // semi-minor axis radius
#define WGS84_1_f       298.257223563  // reciprocal of flattening 1/f
#define WGS84_f         (1.0 / 298.257223563)  // flattening
#define WGS84_e         8.18191908426e-2 // ellipsoidal eccentricity
#define WGS84_e_sq      0.00669437999014 // elliposoidal eccentricity squared
#define WGS84_e12       (WGS84_e_sq / (1 - WGS84_e_sq))
#define WGS84_es        WGS84_e
#define WGS84_ell       (-WGS84_e12)
#define WGS84_Ec        1.5734395859500132
#define WGS84_K2        (-0.0067394967422764341)
#define WGS84_KP2       (1.0067394967422765)
#define WGS84_EPS       (-0.0016792203863837047)
#define WGS84_QUARTER_MERIDIAN (WGS84_b_METERS * WGS84_Ec)

// Max depth required for sncndn. Probably 5 is enough.
#define SNCNDN_DEPTH 5

// common
static inline double eatanhe(double x, double es)
{
    return es > 0.0 ? es * atanh(es * x) : -es * atan(es * x);
}

// common
static inline double taupf(double tau, double es)
{
    double tau1 = hypot(1.0, tau);
    double sig = sinh(eatanhe(tau / tau1, es));
    return hypot(1.0, sig) * tau - sig * tau1;
}

// common
static inline double AngNormalize(double x)
{
    x = fmod(x, 360.0);
    return x < -180 ? x + 360 : (x < 180 ? x : x - 360);
}

// common
static inline double Datan(double x, double y)
{
    double d = x - y;
    double xy = x * y;
    return d ? (2 * xy > -1 ? atan(d / (1 + xy)) : atan(x) - atan(y)) / d :
        1 / (1 + xy);
}

// common
static double SinCosSeries(
    double x, double y, const double c[], int n)
{
    // N.B. n >= 0 and c[] has n+1 elements 0..n, of which c[0] is ignored.
    //
    // Use Clenshaw summation to evaluate
    //   m = (g(x) + g(y)) / 2         -- mean value
    //   s = (g(x) - g(y)) / (x - y)   -- average slope
    // where
    //   g(x) = sum(c[j]*SC(2*j*x), j = 1..n)
    //   SC = sinp ? sin : cos
    //   CS = sinp ? cos : sin
    //
    // This function returns only s; m is discarded.
    //
    // Write
    //   t = [m; s]
    //   t = sum(c[j] * f[j](x,y), j = 1..n)
    // where
    //   f[j](x,y) = [ (SC(2*j*x)+SC(2*j*y))/2 ]
    //               [ (SC(2*j*x)-SC(2*j*y))/d ]
    //
    //             = [       cos(j*d)*SC(j*p)    ]
    //               [ +/-(2/d)*sin(j*d)*CS(j*p) ]
    // (+/- = sinp ? + : -) and
    //    p = x+y, d = x-y
    //
    //   f[j+1](x,y) = A * f[j](x,y) - f[j-1](x,y)
    //
    //   A = [  2*cos(p)*cos(d)      -sin(p)*sin(d)*d]
    //       [ -4*sin(p)*sin(d)/d   2*cos(p)*cos(d)  ]
    //
    // Let b[n+1] = b[n+2] = [0 0; 0 0]
    //     b[j] = A * b[j+1] - b[j+2] + c[j] * I for j = n..1
    //    t =  (c[0] * I  - b[2]) * f[0](x,y) + b[1] * f[1](x,y)
    // c[0] is not accessed for s = t[2]
    double p = x + y, d = x - y,
        cp = cos(p), cd = cos(d),
        sp = sin(p), sd = d ? sin(d) / d : 1,
        m = 2 * cp * cd, s = sp * sd;
    // 2x2 matrices stored in row-major order
    const double a[4] = { m, -s * d * d, -4 * s, m };
    double ba[4] = { 0, 0, 0, 0 };
    double bb[4] = { 0, 0, 0, 0 };
    double* b1 = ba;
    double* b2 = bb;
    if (n > 0) b1[0] = b1[3] = c[n];
    for (int j = n - 1; j > 0; --j) // j = n-1 .. 1
    {
        {
            double* tmp = b1;
            b1 = b2;
            b2 = tmp;
        }
        // b1 = A * b2 - b1 + c[j] * I
        b1[0] = a[0] * b2[0] + a[1] * b2[2] - b1[0] + c[j];
        b1[1] = a[0] * b2[1] + a[1] * b2[3] - b1[1];
        b1[2] = a[2] * b2[0] + a[3] * b2[2] - b1[2];
        b1[3] = a[2] * b2[1] + a[3] * b2[3] - b1[3] + c[j];
    }

    // Here are the full expressions for m and s
    // m =   (c[0] - b2[0]) * f01 - b2[1] * f02 + b1[0] * f11 + b1[1] * f12;
    // s = - b2[2] * f01 + (c[0] - b2[3]) * f02 + b1[2] * f11 + b1[3] * f12;

    // real f01 = 0, f02 = 0;
    double f11 = cd * sp, f12 = 2 * sd * cp;
    // m = b1[0] * f11 + b1[1] * f12;
    s = b1[2] * f11 + b1[3] * f12;

    return s;
}






// direct-only
static inline double max(double a, double b)
{
    return (a > b) ? a : b;
}

// direct-only
static inline double min(double a, double b)
{
    return (a < b) ? a : b;
}

// direct-only
static inline double Delta(double sn, double cn)
{
    return sqrt(WGS84_K2 < 0 ? 1 - WGS84_K2 * sn*sn : WGS84_KP2 + WGS84_K2 * cn*cn);
}

// direct-only
static double RF(double x, double y, double z)
{
    // Carlson, eqs 2.2 - 2.7
    double tolRF = pow(3 * DBL_EPSILON * 0.01, 1.0 / 8.0);
    double
        A0 = (x + y + z) / 3,
        An = A0,
        Q = max(max(fabs(A0 - x), fabs(A0 - y)), fabs(A0 - z)) / tolRF,
        x0 = x,
        y0 = y,
        z0 = z,
        mul = 1;
    while (Q >= mul * fabs(An)) {
        // Max 6 trips
        double lam = sqrt(x0)*sqrt(y0) + sqrt(y0)*sqrt(z0) + sqrt(z0)*sqrt(x0);
        An = (An + lam) / 4;
        x0 = (x0 + lam) / 4;
        y0 = (y0 + lam) / 4;
        z0 = (z0 + lam) / 4;
        mul *= 4;
    }
    double
        X = (A0 - x) / (mul * An),
        Y = (A0 - y) / (mul * An),
        Z = -(X + Y),
        E2 = X*Y - Z*Z,
        E3 = X*Y*Z;
    // http://dlmf.nist.gov/19.36.E1
    // Polynomial is
    // (1 - E2/10 + E3/14 + E2^2/24 - 3*E2*E3/44
    //    - 5*E2^3/208 + 3*E3^2/104 + E2^2*E3/16)
    // convert to Horner form...
    return (E3 * (6930 * E3 + E2 * (15015 * E2 - 16380) + 17160) +
        E2 * ((10010 - 5775 * E2) * E2 - 24024) + 240240) /
        (240240 * sqrt(An));
}

// direct-only
static double RD(double x, double y, double z)
{
    // Carlson, eqs 2.28 - 2.34
    double tolRD = pow(0.2 * (DBL_EPSILON * 0.01), 1.0 / 8.0);
    double
        A0 = (x + y + 3 * z) / 5,
        An = A0,
        Q = max(max(fabs(A0 - x), fabs(A0 - y)), fabs(A0 - z)) / tolRD,
        x0 = x,
        y0 = y,
        z0 = z,
        mul = 1,
        s = 0;
    while (Q >= mul * fabs(An)) {
        // Max 7 trips
        double lam = sqrt(x0)*sqrt(y0) + sqrt(y0)*sqrt(z0) + sqrt(z0)*sqrt(x0);
        s += 1 / (mul * sqrt(z0) * (z0 + lam));
        An = (An + lam) / 4;
        x0 = (x0 + lam) / 4;
        y0 = (y0 + lam) / 4;
        z0 = (z0 + lam) / 4;
        mul *= 4;
    }
    double
        X = (A0 - x) / (mul * An),
        Y = (A0 - y) / (mul * An),
        Z = -(X + Y) / 3,
        E2 = X*Y - 6 * Z*Z,
        E3 = (3 * X*Y - 8 * Z*Z)*Z,
        E4 = 3 * (X*Y - Z*Z) * Z*Z,
        E5 = X*Y*Z*Z*Z;
    // http://dlmf.nist.gov/19.36.E2
    // Polynomial is
    // (1 - 3*E2/14 + E3/6 + 9*E2^2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26
    //    - E2^3/16 + 3*E3^2/40 + 3*E2*E4/20 + 45*E2^2*E3/272
    //    - 9*(E3*E4+E2*E5)/68)
    return ((471240 - 540540 * E2) * E5 +
        (612612 * E2 - 540540 * E3 - 556920) * E4 +
        E3 * (306306 * E3 + E2 * (675675 * E2 - 706860) + 680680) +
        E2 * ((417690 - 255255 * E2) * E2 - 875160) + 4084080) /
        (4084080 * mul * An * sqrt(An)) + 3 * s;
}

// direct-only
static double E(double sn, double cn, double dn)
{
    double
        cn2 = cn*cn, dn2 = dn*dn, sn2 = sn*sn,
        ei = RF(cn2, dn2, 1) - WGS84_K2 * sn2 * RD(cn2, dn2, 1) / 3;
    //ei = (_k2 <= 0 ?
    //    // Carlson, eq. 4.6 and
    //    // http://dlmf.nist.gov/19.25.E9
    //    RF(cn2, dn2, 1) - _k2 * sn2 * RD(cn2, dn2, 1) / 3 :
    //    (_kp2 >= 0 ?
    //        // http://dlmf.nist.gov/19.25.E10
    //        _kp2 * RF(cn2, dn2, 1) +
    //        _k2 * _kp2 * sn2 * RD(cn2, 1, dn2) / 3 +
    //        _k2 * abs(cn) / dn :
    //        // http://dlmf.nist.gov/19.25.E11
    //        -_kp2 * sn2 * RD(dn2, 1, cn2) / 3 + dn / abs(cn)));
    ei *= fabs(sn);
    // Enforce usual trig-like symmetries
    if (cn < 0)
        ei = 2 * WGS84_Ec - ei;
    if (sn < 0)
        ei = -ei;
    return ei;
}

// direct-only
static inline double Ed(double ang)
{
    double sn = sin(ang);
    double cn = cos(ang);
    return E(sn, cn, Delta(sn, cn));
}

// direct-only
static double Einv(double x)
{
    double tolJAC = sqrt(DBL_EPSILON * 0.01);
    double n = floor(x / (2 * WGS84_Ec) + 0.5);
    x -= 2 * WGS84_Ec * n;      // x now in [-ec, ec)
                                // Linear approximation
    double phi = PI * x / (2 * WGS84_Ec); // phi in [-pi/2, pi/2)
                                          // First order correction
    phi -= WGS84_EPS * sin(2 * phi) / 2;
    for (int i = 0; i < SNCNDN_DEPTH; ++i) {
        double
            sn = sin(phi),
            cn = cos(phi),
            dn = Delta(sn, cn),
            err = (E(sn, cn, dn) - x) / dn;
        phi = phi - err;
        if (fabs(err) < tolJAC)
            break;
    }
    return n * PI + phi;
}

// direct-only
static inline double InverseRectifyingLatitude(double mu)
{
    if (fabs(mu) == (PI / 2.0))
        return mu;

    double beta = Einv(mu * WGS84_Ec / (PI / 2.0));
    double inverse_parametric_lat = atan(tan(beta) / (1.0 - WGS84_f));
    return inverse_parametric_lat;
}

// direct-only
static inline double Dasinh(double x, double y)
{
    double d = x - y;
    double hx = hypot(1.0, x);
    double hy = hypot(1.0, y);
    return d ? asinh(x*y > 0 ? d * (x + y) / (x*hy + y*hx) :
        x*hy - y*hx) / d :
        1 / hx;
}

// direct-only
static inline double Dgdinv(double x, double y)
{
    return Dasinh(x, y) / Datan(x, y);
}

// direct-only
static inline double DRectifyingToConformal(double mux, double muy)
{
    const int TM_MAXORD = 6; // from a precision of 2 in GeographicLib

    static const double bet[] = {
        0.00000000000000000,
        0.00083773216405794864,
        5.9058701522202026e-08,
        1.6734826652839968e-10,
        2.1647980400627056e-13,
        3.7879780461686058e-16,
        7.2487488906941545e-19,
    };

    return 1 - SinCosSeries(mux, muy, bet, TM_MAXORD);
}


// direct-only
static inline double DRectifyingToIsometric(double mux, double muy)
{
    double latx = InverseRectifyingLatitude(mux);
    double laty = InverseRectifyingLatitude(muy);
    return
        Dgdinv(taupf(tan(latx), WGS84_es),
               taupf(tan(laty), WGS84_es)) *
        DRectifyingToConformal(mux, muy);
}

// direct-only
RhumbLine RhumbLine_create(double lat1, double lon1, double heading)
{
    RhumbLine line;

    // Fudge the starting point at the poles.
    lat1 = min(lat1, 89.999999);
    lat1 = max(lat1, -89.999999);

    double lat1_rad = DEG_TO_RAD(lat1);

    line.lat1 = lat1;
    line.lon1 = lon1;
    line.azi12 = AngNormalize(heading); // [-180, 180)

    double alp12 = DEG_TO_RAD(line.azi12);
    line.salp = line.azi12 == -180 ? 0 : sin(alp12);
    line.calp = fabs(line.azi12) == 90 ? 0 : cos(alp12);

    // get the rectifying latitude
    double parametric_lat = atan((1.0 - WGS84_f) * tan(lat1_rad));
    double meridian_dist = WGS84_b_METERS * Ed(parametric_lat);
    line.mu1 = 90 * meridian_dist / WGS84_QUARTER_MERIDIAN;

    // get the circle radius at the specified latitude
    double f1 = 1.0 - (1.0f / WGS84_1_f);
    // a * cos(beta)
    line.r1 = WGS84_a_METERS / hypot(1.0, f1 * tan(lat1_rad));

    return line;
}

// direct-only
void RhumbLine_Direct(
    double lat1, double lon1, double heading, double dist,
    double* lat2, double* lon2)
{
    RhumbLine line = RhumbLine_create(lat1, lon1, heading);

    double mu12 = dist * line.calp * 90 / WGS84_QUARTER_MERIDIAN;
    double mu2 = line.mu1 + mu12;
    double lat2x, lon2x;
    if (fabs(mu2) <= 90)
    {
        if (line.calp)
        {
            lat2x = RAD_TO_DEG(InverseRectifyingLatitude(DEG_TO_RAD(mu2)));
            double psi12 = DRectifyingToIsometric(
                DEG_TO_RAD(mu2), DEG_TO_RAD(line.mu1)) * mu12;
            lon2x = line.salp * psi12 / line.calp;
        }
        else
        {
            lat2x = line.lat1;
            lon2x = line.salp * dist / line.r1;
        }
        lon2x = AngNormalize(AngNormalize(line.lon1) + lon2x);
    }
    else
    {
        // Reduce to the interval [-180, 180)
        mu2 = AngNormalize(mu2);
        // Deal with points on the anti-meridian
        if (fabs(mu2) > 90) mu2 = AngNormalize(180 - mu2);
        lat2x = RAD_TO_DEG(InverseRectifyingLatitude(DEG_TO_RAD(mu2)));
        lon2x = NAN;
    }

    *lat2 = lat2x;
    *lon2 = lon2x;
}


















// inverse-only
double DConformalToRectifying(double chix, double chiy)
{
    const int TM_MAXORD = 6; // from a precision of 2 in GeographicLib

    static const double alp[] = {
        0.00000000000000000,
        0.00083773182062446983,
        7.6085277735723085e-07,
        1.1976455033294525e-09,
        2.4291706072013591e-12,
        5.7117576778658038e-15,
        1.4911177312583895e-17,
    };

    return 1 + SinCosSeries(chix, chiy, alp, TM_MAXORD);
}

// inverse-only
static inline double gd(double x)
{
    return atan(sinh(x));
}

// inverse-only
static inline double Dsinh(double x, double y) {
    double d = (x - y) / 2;
    return cosh((x + y) / 2) * (d ? sinh(d) / d : 1);
}

// inverse-only
static inline double Dgd(double x, double y) {
    return Datan(sinh(x), sinh(y)) * Dsinh(x, y);
}

// inverse-only
static inline double sum(double u, double v, double* t)
{
    volatile double s = u + v;
    volatile double up = s - v;
    volatile double vpp = s - up;
    up -= u;
    vpp -= v;
    *t = -(up + vpp);
    // u + v =       s      + t
    //       = round(u + v) + t
    return s;
}

// inverse-only
static inline double AngDiff(double x, double y)
{
    double t;
    double d = -AngNormalize(
        sum(remainder(x, 360.0), remainder(-y, 360.0), &t));

    // Here y - x = d - t (mod 360), exactly, where d is in (-180,180] and
    // abs(t) <= eps (eps = 2^-45 for doubles).  The only case where the
    // addition of t takes the result outside the range (-180,180] is d = 180
    // and t < 0.  The case, d = -180 + eps, t = eps, can't happen, since
    // sum would have returned the exact result in such a case (i.e., given t
    // = 0).
    return (d == 180 && t < 0 ? -180 : d) - t;
}

// inverse-only
static inline double IsometricLatitude(double phi)
{
    phi = DEG_TO_RAD(phi);
    return RAD_TO_DEG(asinh(taupf(tan(phi), WGS84_es)));
}

// inverse-only
double DIsometricToRectifying(double psix, double psiy)
{
    psix = DEG_TO_RAD(psix);
    psiy = DEG_TO_RAD(psiy);
    return DConformalToRectifying(gd(psix), gd(psiy)) * Dgd(psix, psiy);
}

// inverse-only
void RhumbLine_Inverse(
    double lat1, double lon1, double lat2, double lon2,
    double* dist, double* heading)
{
    double
        lon12 = AngDiff(lon1, lon2),
        psi1 = IsometricLatitude(lat1),
        psi2 = IsometricLatitude(lat2),
        psi12 = psi2 - psi1,
        h = hypot(lon12, psi12);
    *heading = RAD_TO_DEG(atan2(lon12, psi12));
    double dmudpsi = DIsometricToRectifying(psi2, psi1);
    *dist = h * dmudpsi * WGS84_QUARTER_MERIDIAN / 90;
}

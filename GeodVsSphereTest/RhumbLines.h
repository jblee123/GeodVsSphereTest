#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    double lat1, lon1;
    double azi12;  // normalized heading in degrees [-180, 180)
    double salp;  // sin(heading)
    double calp;  // cos(heading)
    double mu1;  // rectifying lat (in degrees)
    double r1;  // circle radius at lat1
} RhumbLine;

void RhumbLine_Direct(
    double lat1, double lon1, double heading, double dist,
    double* lat2, double* lon2);

void RhumbLine_Inverse(
    double lat1, double lon1, double lat2, double lon2,
    double* dist, double* heading);

#ifdef __cplusplus
}
#endif

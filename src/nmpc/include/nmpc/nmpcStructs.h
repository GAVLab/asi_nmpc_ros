#ifndef NMPCSTRUCTS_H
#define NMPCSTRUCTS_H
#include <vector>

struct nmpcFlags
{
    bool newGpsPos;
    bool newGpsCourse;
    bool newVel;
    bool newYawRate;

    bool init;
    bool reinit;

};

struct gpsData
{
    double lat;
    double lon;
    double alt;
    double course;
    double velMag;
    double localVel_long;
    std::vector <double> nedVel;
};

struct imuData
{
    double yawRate;
};

#endif // NMPCSTRUCTS_H

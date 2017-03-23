
#include "GpsFunctions.h"

double GpsFunctions::checkWrap(double angle) //double
{
//    while (angle<-PI) {
//        angle += (2. * PI);
//    }
//
//    while (angle >= PI) {
//        angle -= (2. * PI);
//    }

    return (fmod(angle+PI,2*PI)-PI);
}

int GpsFunctions::timeStamp(int &gps_weeks, int &gps_ms) //int
{
    //int initial_time;
    //int present_time;
    //int last_time;
    int wk_to_ms;
    int num_weeks;
    int ms;
    ptime initial_time(date(1980, Jan, 6));
    ptime present_time(microsec_clock::universal_time());
    time_duration diff = present_time - initial_time;
    wk_to_ms = 7 * 24 * 60 * 60 * 1E3;
    num_weeks = (diff.total_milliseconds()) / wk_to_ms;
    ms = (diff.total_milliseconds() % wk_to_ms) + 15000; //add 15s for gps time
    gps_weeks = num_weeks;
    gps_ms = ms;

    return (EXIT_SUCCESS);
}
//END GpsFunciton namespace
double GpsFunctions::check_t(double time) {
    //CHECK_T accounting for beginning or end of week crossover.
    const double half_week = 302400; //seconds
    double corrTime = time;
    if (time > half_week)
        corrTime = time - 2 * half_week;
    else if (time < -half_week)
        corrTime = time + 2 * half_week;
    return corrTime;

} //end check_t

double GpsFunctions::wrap_heading(double heading) {
    double x = heading;
    while (x >= gpsPi) {
        x = x - 2 * gpsPi;
    }
    while (x<-gpsPi) {
        x = x + 2 * gpsPi;
    }
    return (x);
}

int GpsFunctions::wgslla2utm(double Lat, double Long, double &northing, double &easting, int &zone, bool &north) {
    const double a = 6378137.0;
    const double ee = 0.00669437999;
    const double k0 = 0.9996;
    const double e2 = ee / (1 - ee);

    double LongTemp = (Long + 180) - int((Long + 180) / 360)*360 - 180; // -180.00 .. 179.9;
    double LatRad = GRAD_A_RAD(Lat);
    double LongRad = GRAD_A_RAD(LongTemp);
    double LongOriginRad;

    double N, T, C, A, M;

    //Make sure the longitude is between -180.00 .. 179.9
    zone = int((LongTemp + 180) / 6.0) + 1;
    if (Lat >= 56.0 && Lat < 64.0 && LongTemp >= 3.0 && LongTemp < 12.0)
        zone = 32;

    // Special zones for Svalbard
    if (Lat >= 72.0 && Lat < 84.0) {
        if (LongTemp >= 0.0 && LongTemp < 9.0)
            zone = 31;
        else if (LongTemp >= 9.0 && LongTemp < 21.0)
            zone = 33;
        else if (LongTemp >= 21.0 && LongTemp < 33.0)
            zone = 35;
        else if (LongTemp >= 33.0 && LongTemp < 42.0)
            zone = 37;
    }
    LongOriginRad = GRAD_A_RAD((zone - 1)*6 - 180 + 3);

    N = a / sqrt(1 - ee * sin(LatRad) * sin(LatRad));
    T = tan(LatRad) * tan(LatRad);
    C = e2 * cos(LatRad) * cos(LatRad);
    A = cos(LatRad)*(LongRad - LongOriginRad);
    M = a * ((1 - ee / 4 - 3 * ee * ee / 64 - 5 * ee * ee * ee / 256) * LatRad
            - (3 * ee / 8 + 3 * ee * ee / 32 + 45 * ee * ee * ee / 1024) * sin(2 * LatRad)
            + (15 * ee * ee / 256 + 45 * ee * ee * ee / 1024) * sin(4 * LatRad)
            - (35 * ee * ee * ee / 3072) * sin(6 * LatRad));

    easting = (double) (k0 * N * (A + (1 - T + C) * A * A * A / 6
            + (5 - 18 * T + T * T + 72 * C - 58 * e2) * A * A * A * A * A / 120) + 500000.0);
    northing = (double) (k0 * (M + N * tan(LatRad)*(A * A / 2 + (5 - T + 9 * C + 4 * C * C) * A * A * A * A / 24
            + (61 - 58 * T + T * T + 600 * C - 330 * e2) * A * A * A * A * A * A / 720)));

    if (Lat < 0) {
        northing += 10000000; //10000000 meter offset for southern hemisphere
        north = false;
    } else
        north = true;

    return 1;
}

int GpsFunctions::wgslla2xyz(double xyz[], double wlat, double wlon, double walt) {
    double A_EARTH = 6378137.;
    double flattening = 1. / 298.257223563;
    double NAV_E2 = (2.0 - flattening) * flattening; // also e^2
    double deg2rad = gpsPi / 180.;

    double slat = sin(wlat * deg2rad);
    double clat = cos(wlat * deg2rad);
    double r_n = A_EARTH / sqrt(1 - NAV_E2 * slat * slat);
    xyz[0] = (r_n + walt) * clat * cos(wlon * deg2rad);
    xyz[1] = (r_n + walt) * clat * sin(wlon * deg2rad);
    xyz[2] = (r_n * (1 - NAV_E2) + walt) * slat;

    if ((wlat < -90.0) | (wlat > +90.0) | (wlon < -180.0) | (wlon > +360.0)) {
        cout << "GpsFunctions.cc:wgslla2xyz-> WGS lat or WGS lon out of range in wgslla2xyz.cpp\r" << endl;
        return (0);
    } else {
        return (1);
    }
}

int GpsFunctions::rot(double R[3][3], double angle, int axis) {

    double cang = cos(angle * gpsPi / 180);
    double sang = sin(angle * gpsPi / 180);

    if (axis == 1) {
        R[0][0] = 1;
        R[0][1] = 0;
        R[0][2] = 0;
        R[1][0] = 0;
        R[2][0] = 0;
        R[1][1] = cang;
        R[2][2] = cang;
        R[1][2] = sang;
        R[2][1] = -sang;
    } else if (axis == 2) {
        R[0][1] = 0;
        R[1][0] = 0;
        R[1][1] = 1;
        R[1][2] = 0;
        R[2][1] = 0;
        R[0][0] = cang;
        R[2][2] = cang;
        R[0][2] = -sang;
        R[2][0] = sang;
    } else if (axis == 3) {
        R[2][0] = 0;
        R[2][1] = 0;
        R[2][2] = 1;
        R[0][2] = 0;
        R[1][2] = 0;
        R[0][0] = cang;
        R[1][1] = cang;
        R[1][0] = -sang;
        R[0][1] = sang;
    }
    return (1);
}

int GpsFunctions::rot(ublas::matrix<double> &R, double angle, int axis) {

    double cang = cos(angle * gpsPi / 180);
    double sang = sin(angle * gpsPi / 180);

    if (axis == 1) {
        R(0,0) = 1;
        R(0,1) = 0;
        R(0,2) = 0;
        R(1,0) = 0;
        R(2,0) = 0;
        R(1,1) = cang;
        R(2,2) = cang;
        R(1,2) = sang;
        R(2,1) = -sang;
    } else if (axis == 2) {
        R(0,1) = 0;
        R(1,0) = 0;
        R(1,1) = 1;
        R(1,2) = 0;
        R(2,1) = 0;
        R(0,0) = cang;
        R(2,2) = cang;
        R(0,2) = -sang;
        R(2,0) = sang;
    } else if (axis == 3) {
        R(2,0) = 0;
        R(2,1) = 0;
        R(2,2) = 1;
        R(0,2) = 0;
        R(1,2) = 0;
        R(0,0) = cang;
        R(1,1) = cang;
        R(1,0) = -sang;
        R(0,1) = sang;
    }
    return (1);
}

int GpsFunctions::wgsxyz2enu(double enu[], double xyz[], double reflat, double reflon, double refalt) {
    //First, calculate the xyz of reflat, reflon, refalt
    double refxyz[3];
    int lla2xyz_ret = wgslla2xyz(refxyz, reflat, reflon, refalt);

    //Difference xyz from reference point
    double diffxyz[3];
    diffxyz[0] = xyz[0] - refxyz[0];
    diffxyz[1] = xyz[1] - refxyz[1];
    diffxyz[2] = xyz[2] - refxyz[2];

    // Now rotate the (often short) diffxyz vector to enu frame
    double R1[3][3];
    double R2[3][3];
    int rot_ret1 = rot(R1, 90 + reflon, 3);
    int rot_ret2 = rot(R2, 90 - reflat, 1);
    double R[3][3];
    R[0][0] = R2[0][0] * R1[0][0] + R2[0][1] * R1[1][0] + R2[0][2] * R1[2][0];
    R[0][1] = R2[0][0] * R1[0][1] + R2[0][1] * R1[1][1] + R2[0][2] * R1[2][1];
    R[0][2] = R2[0][0] * R1[0][2] + R2[0][1] * R1[1][2] + R2[0][2] * R1[2][2];
    R[1][0] = R2[1][0] * R1[0][0] + R2[1][1] * R1[1][0] + R2[1][2] * R1[2][0];
    R[1][1] = R2[1][0] * R1[0][1] + R2[1][1] * R1[1][1] + R2[1][2] * R1[2][1];
    R[1][2] = R2[1][0] * R1[0][2] + R2[1][1] * R1[1][2] + R2[1][2] * R1[2][2];
    R[2][0] = R2[2][0] * R1[0][0] + R2[2][1] * R1[1][0] + R2[2][2] * R1[2][0];
    R[2][1] = R2[2][0] * R1[0][1] + R2[2][1] * R1[1][1] + R2[2][2] * R1[2][1];
    R[2][2] = R2[2][0] * R1[0][2] + R2[2][1] * R1[1][2] + R2[2][2] * R1[2][2];

    enu[0] = R[0][0] * diffxyz[0] + R[0][1] * diffxyz[1] + R[0][2] * diffxyz[2];
    enu[1] = R[1][0] * diffxyz[0] + R[1][1] * diffxyz[1] + R[1][2] * diffxyz[2];
    enu[2] = R[2][0] * diffxyz[0] + R[2][1] * diffxyz[1] + R[2][2] * diffxyz[2];

    //cout << "R[0][0]= " << R[0][0] << "\tR[0][1]= " << R[0][1] << "\tR[0][2]= " << R[0][2] << "\r" << endl;

    return (1);
}

int GpsFunctions::wgsdiffxyz2diffenu(double diffenu[], double diffxyz[], double reflat, double reflon) {
    //First, calculate the xyz of reflat, reflon, refalt
    //        double refxyz[3];
    //        int lla2xyz_ret = wgslla2xyz(refxyz, reflat, reflon, refalt);

    //Difference xyz from reference point
    //        double diffxyz[3];
    //        diffxyz[0] = xyz[0] - refxyz[0];
    //        diffxyz[1] = xyz[1] - refxyz[1];
    //        diffxyz[2] = xyz[2] - refxyz[2];

    // Now rotate the (often short) diffxyz vector to enu frame
    double R1[3][3];
    double R2[3][3];
    int rot_ret1 = rot(R1, 90 + reflon, 3);
    int rot_ret2 = rot(R2, 90 - reflat, 1);
    double R[3][3];

    //[ R2_00*R1_00+R2_01*R1_10+R2_02*R1_20, R2_00*R1_01+R2_01*R1_11+R2_02*R1_21, R2_00*R1_02+R2_01*R1_12+R2_02*R1_22]
    //[ R2_10*R1_00+R2_11*R1_10+R2_12*R1_20, R2_10*R1_01+R2_11*R1_11+R2_12*R1_21, R2_10*R1_02+R2_11*R1_12+R2_12*R1_22]
    //[ R2_20*R1_00+R2_21*R1_10+R2_22*R1_20, R2_20*R1_01+R2_21*R1_11+R2_22*R1_21, R2_20*R1_02+R2_21*R1_12+R2_22*R1_22]


    R[0][0] = R2[0][0] * R1[0][0] + R2[0][1] * R1[1][0] + R2[0][2] * R1[2][0];
    R[1][0] = R2[1][0] * R1[0][0] + R2[1][1] * R1[1][0] + R2[1][2] * R1[2][0];
    R[2][0] = R2[2][0] * R1[0][0] + R2[2][1] * R1[1][0] + R2[2][2] * R1[2][0];

    R[0][1] = R2[0][0] * R1[0][1] + R2[0][1] * R1[1][1] + R2[0][2] * R1[2][1];
    R[1][1] = R2[1][0] * R1[0][1] + R2[1][1] * R1[1][1] + R2[1][2] * R1[2][1];
    R[2][1] = R2[2][0] * R1[0][1] + R2[2][1] * R1[1][1] + R2[2][2] * R1[2][1];

    R[0][2] = R2[0][0] * R1[0][2] + R2[0][1] * R1[1][2] + R2[0][2] * R1[2][2];
    R[1][2] = R2[1][0] * R1[0][2] + R2[1][1] * R1[1][2] + R2[1][2] * R1[2][2];
    R[2][2] = R2[2][0] * R1[0][2] + R2[2][1] * R1[1][2] + R2[2][2] * R1[2][2];

    diffenu[0] = R[0][0] * diffxyz[0] + R[0][1] * diffxyz[1] + R[0][2] * diffxyz[2];
    diffenu[1] = R[1][0] * diffxyz[0] + R[1][1] * diffxyz[1] + R[1][2] * diffxyz[2];
    diffenu[2] = R[2][0] * diffxyz[0] + R[2][1] * diffxyz[1] + R[2][2] * diffxyz[2];

    //cout << "R[0][0]= " << R[0][0] << "\tR[0][1]= " << R[0][1] << "\tR[0][2]= " << R[0][2] << "\r" << endl;

    return (1);
}

int GpsFunctions::wgsdiffxyz2diffenu(ublas::vector<double> &diffenu, std::vector<double> &diffxyz, double reflat, double reflon) {
    //First, calculate the xyz of reflat, reflon, refalt
    //        double refxyz[3];
    //        int lla2xyz_ret = wgslla2xyz(refxyz, reflat, reflon, refalt);

    //Difference xyz from reference point
    //        double diffxyz[3];
    //        diffxyz[0] = xyz[0] - refxyz[0];
    //        diffxyz[1] = xyz[1] - refxyz[1];
    //        diffxyz[2] = xyz[2] - refxyz[2];

    // Now rotate the (often short) diffxyz vector to enu frame

    double R1[3][3];
    double R2[3][3];
    int rot_ret1 = rot(R1, 90 + reflon, 3);
    int rot_ret2 = rot(R2, 90 - reflat, 1);
    double R[3][3];

    //[ R2_00*R1_00+R2_01*R1_10+R2_02*R1_20, R2_00*R1_01+R2_01*R1_11+R2_02*R1_21, R2_00*R1_02+R2_01*R1_12+R2_02*R1_22]
    //[ R2_10*R1_00+R2_11*R1_10+R2_12*R1_20, R2_10*R1_01+R2_11*R1_11+R2_12*R1_21, R2_10*R1_02+R2_11*R1_12+R2_12*R1_22]
    //[ R2_20*R1_00+R2_21*R1_10+R2_22*R1_20, R2_20*R1_01+R2_21*R1_11+R2_22*R1_21, R2_20*R1_02+R2_21*R1_12+R2_22*R1_22]


    R[0][0] = R2[0][0] * R1[0][0] + R2[0][1] * R1[1][0] + R2[0][2] * R1[2][0];
    R[1][0] = R2[1][0] * R1[0][0] + R2[1][1] * R1[1][0] + R2[1][2] * R1[2][0];
    R[2][0] = R2[2][0] * R1[0][0] + R2[2][1] * R1[1][0] + R2[2][2] * R1[2][0];

    R[0][1] = R2[0][0] * R1[0][1] + R2[0][1] * R1[1][1] + R2[0][2] * R1[2][1];
    R[1][1] = R2[1][0] * R1[0][1] + R2[1][1] * R1[1][1] + R2[1][2] * R1[2][1];
    R[2][1] = R2[2][0] * R1[0][1] + R2[2][1] * R1[1][1] + R2[2][2] * R1[2][1];

    R[0][2] = R2[0][0] * R1[0][2] + R2[0][1] * R1[1][2] + R2[0][2] * R1[2][2];
    R[1][2] = R2[1][0] * R1[0][2] + R2[1][1] * R1[1][2] + R2[1][2] * R1[2][2];
    R[2][2] = R2[2][0] * R1[0][2] + R2[2][1] * R1[1][2] + R2[2][2] * R1[2][2];

    diffenu(0) = R[0][0] * diffxyz[0] + R[0][1] * diffxyz[1] + R[0][2] * diffxyz[2];
    diffenu(1) = R[1][0] * diffxyz[0] + R[1][1] * diffxyz[1] + R[1][2] * diffxyz[2];
    diffenu(2) = R[2][0] * diffxyz[0] + R[2][1] * diffxyz[1] + R[2][2] * diffxyz[2];

    //cout << "R[0][0]= " << R[0][0] << "\tR[0][1]= " << R[0][1] << "\tR[0][2]= " << R[0][2] << "\r" << endl;

    return (1);
}

int GpsFunctions::wgslla2enu(double enu[], double lat, double lon, double alt, double reflat, double reflon, double refalt) {
    double xyz[3];
    int lla2xyz_ret = wgslla2xyz(xyz, lat, lon, alt);
    int xyz2enu_ret = wgsxyz2enu(enu, xyz, reflat, reflon, refalt);

    return (1);
}
/*
int GpsFunctions::wgsxyz2lla(double &lat, double &lon, double &alt, double xyz[]) {

    //This dual-variable iteration seems to be 7 or 8 times faster than
    //a one-variable (in latitude only) iteration.  AKB 7/17/95

    double A_EARTH = 6378137;
    double flattening = 1 / 298.257223563;
    double NAV_E2 = (2 - flattening) * flattening; // also e^2
    double rad2deg = 180 / gpsPi;

    if ((xyz[0] == 0.0) & (xyz[1] == 0.0)) {
        lon = 0.0;
    } else {
        lon = atan2(xyz[1], xyz[0]) * rad2deg;
    }

    if ((xyz[0] == 0.0) & (xyz[1] == 0.0) & (xyz[2] == 0.0)) {
        //error('WGS xyz at center of earth');

    } else {
        // Make initial lat and alt guesses based on spherical earth.
        double rhosqrd = xyz[0] * xyz[0] + xyz[1] * xyz[1];
        double rho = sqrt(rhosqrd);
        double templat = atan2(xyz[2], rho);
        double tempalt = sqrt(rhosqrd + xyz[2] * xyz[2]) - A_EARTH;
        double rhoerror = 1000.0;
        double zerror = 1000.0;

        //		%  Newton's method iteration on templat and tempalt makes
        //		%	the residuals on rho and z progressively smaller.  Loop
        //		%	is implemented as a 'while' instead of a 'do' to simplify
        //		%	porting to MATLAB

        while ((abs(rhoerror) > 1e-6) | (abs(zerror) > 1e-6)) {
            double slat = sin(templat);
            double clat = cos(templat);
            double q = 1 - NAV_E2 * slat*slat;
            double r_n = A_EARTH / sqrt(q);
            double drdl = r_n * NAV_E2 * slat * clat / q; // d(r_n)/d(latitutde)

            rhoerror = (r_n + tempalt) * clat - rho;
            zerror = (r_n * (1 - NAV_E2) + tempalt) * slat - xyz[2];

            //			%			  --                               -- --      --
            //			%			  |  drhoerror/dlat  drhoerror/dalt | |  a  b  |
            //                        % Find Jacobian           |		       		    |=|        |
            //			%			  |   dzerror/dlat    dzerror/dalt  | |  c  d  |
            //			%			  --                               -- --      --

            double aa = drdl * clat - (r_n + tempalt) * slat;
            double bb = clat;
            double cc = (1 - NAV_E2)*(drdl * slat + r_n * clat);
            double dd = slat;

            //Apply correction = inv(Jacobian)*errorvector

            double invdet = 1.0 / (aa * dd - bb * cc);
            templat = templat - invdet * (+dd * rhoerror - bb * zerror);
            tempalt = tempalt - invdet * (-cc * rhoerror + aa * zerror);
        }

        lat = templat*rad2deg;
        alt = tempalt;
    }
    return (1);
}*/

int GpsFunctions::wgsxyz2lla(double &lat, double &lon, double &alt, std::vector<double> &xyz) {

    //This dual-variable iteration seems to be 7 or 8 times faster than
    //a one-variable (in latitude only) iteration.  AKB 7/17/95

    //cout << "Called correct function.\n" << endl;

    double A_EARTH = 6378137;
    double flattening = 1 / 298.257223563;
    double NAV_E2 = (2 - flattening) * flattening; // also e^2
    double rad2deg = 180 / gpsPi;

    //cout << "Righting xyz to the screen.\n" << endl;
    //cout << "Size of xyz.\n" << xyz.size() << endl;
    //cout << "xyz X = " << xyz[0] << endl; 
    if ((xyz[0] == 0.0) & (xyz[1] == 0.0)) {
        lon = 0.0;
    } else {
        lon = atan2(xyz[1], xyz[0]) * rad2deg;
    }

    //cout << "Longitude Assigned.\n" << endl;
    if ((xyz[0] == 0.0) & (xyz[1] == 0.0) & (xyz[2] == 0.0)) {
        cout << "WGS xyz at center of earth.\n" << endl;

    } else {
        // Make initial lat and alt guesses based on spherical earth.
        double rhosqrd = xyz[0] * xyz[0] + xyz[1] * xyz[1];
        double rho = sqrt(rhosqrd);
        double templat = atan2(xyz[2], rho);
        double tempalt = sqrt(rhosqrd + xyz[2] * xyz[2]) - A_EARTH;
        double rhoerror = 1000.0;
        double zerror = 1000.0;

        //      %  Newton's method iteration on templat and tempalt makes
        //      %   the residuals on rho and z progressively smaller.  Loop
        //      %   is implemented as a 'while' instead of a 'do' to simplify
        //      %   porting to MATLAB

        while ((abs(rhoerror) > 1e-6) | (abs(zerror) > 1e-6)) {
            double slat = sin(templat);
            double clat = cos(templat);
            double q = 1 - NAV_E2 * slat*slat;
            double r_n = A_EARTH / sqrt(q);
            double drdl = r_n * NAV_E2 * slat * clat / q; // d(r_n)/d(latitutde)

            rhoerror = (r_n + tempalt) * clat - rho;
            zerror = (r_n * (1 - NAV_E2) + tempalt) * slat - xyz[2];

            //          %             --                               -- --      --
            //          %             |  drhoerror/dlat  drhoerror/dalt | |  a  b  |
            //                        % Find Jacobian           |                       |=|        |
            //          %             |   dzerror/dlat    dzerror/dalt  | |  c  d  |
            //          %             --                               -- --      --

            double aa = drdl * clat - (r_n + tempalt) * slat;
            double bb = clat;
            double cc = (1 - NAV_E2)*(drdl * slat + r_n * clat);
            double dd = slat;

            //Apply correction = inv(Jacobian)*errorvector

            double invdet = 1.0 / (aa * dd - bb * cc);
            templat = templat - invdet * (+dd * rhoerror - bb * zerror);
            tempalt = tempalt - invdet * (-cc * rhoerror + aa * zerror);
        }

        lat = templat*rad2deg;
        alt = tempalt;
    }
    return (1);
}


//get rotation matrix from ECEF to NED
ublas::matrix<double> GpsFunctions::get_C_xyz2ned(ublas::matrix<double> x,double error)
{
	ublas::matrix<double> C(3,3);

	double p = sqrt(x(0,0)*x(0,0)+x(1,0)*x(1,0));

	//initialize Rn and h
	double Rn = eRad;
	double h = 0;
	double e_squared = ecc*ecc;
	//e_squared = .006694379990141316;

	double lon = atan2(x(1,0),x(0,0));

	double min_error = 1;
	double ii = 1;
	double lastRn = 0;
	double lat;

	while (min_error > error)
	{
                double sin_latitude=x(2,0)/((1-e_squared)*Rn+h);
                lat = atan((x(2,0)+e_squared*Rn*sin_latitude)/p);

                lastRn=Rn;
                Rn=eRad/(sqrt((1-e_squared*(sin(lat))*(sin(lat)))));
                min_error=abs(lastRn-Rn);
                h=(p/cos(lat))-Rn;
	}

	C(0,0)=-sin(lat)*cos(lon);
	C(0,1) = -sin(lat)*sin(lon);
	C(0,2) = cos(lat);
	C(1,0) = -sin(lon);
	C(1,1) = cos(lon);
	C(1,2) = 0;
	C(2,0) = -cos(lat)*cos(lon);
	C(2,1) = -cos(lat)*sin(lon);
	C(2,2) = -sin(lat);	

	return C;


}

ublas::matrix<double> GpsFunctions::get_C_wgslla2ned(ublas::matrix<double> x)
{
	ublas::matrix<double> C(3,3);

	double lat = x(0,0);
	double lon = x(1,0);

	C(0,0)=-sin(lat)*cos(lon);
	C(0,1) = -sin(lat)*sin(lon);
	C(0,2) = cos(lat);
	C(1,0) = -sin(lon);
	C(1,1) = cos(lon);
	C(1,2) = 0;
	C(2,0) = -cos(lat)*cos(lon);
	C(2,1) = -cos(lat)*sin(lon);
	C(2,2) = -sin(lat);

	return C;
}

ublas::matrix<double> GpsFunctions::xyz2wgslla(ublas::matrix<double> x,double error)
{
	ublas::matrix<double> xlla(3,1);

	double p = sqrt(x(0,0)*x(0,0)+x(1,0)*x(1,0));

	//initialize Rn and h
	double Rn = eRad;
	double h = 0;
	double e_squared = ecc*ecc;

	double lon = atan2(x(1,0),x(0,0));

	double min_error = 1;
	double ii = 1;
	double lastRn = 0;
	double lat;

	while (min_error > error)
	{
                double sin_latitude=x(2,0)/((1-e_squared)*Rn+h);
                lat = atan((x(2,0)+e_squared*Rn*sin_latitude)/p);

                lastRn=Rn;
                Rn=eRad/(sqrt((1-e_squared*(sin(lat))*(sin(lat)))));
                min_error=abs(lastRn-Rn);
                h=(p/cos(lat))-Rn;
	}

	xlla(0,0) = lat;
	xlla(1,0) = lon;
	xlla(2,0) = h;

	return xlla;
}

ublas::matrix<double> GpsFunctions::wgslla2xyz(ublas::matrix<double> x)
{
	//	% constants
//a=6378137;
//e_squared=6.694379990141316e-003;
//
//%computation
//Rn=a/((1-e_squared*(sin(x(1,1)))^2)^(1/2));
//out(1,1)=(Rn+x(3,1))*cos(x(1,1))*cos(x(2,1));
//out(2,1)=(Rn+x(3,1))*cos(x(1,1))*sin(x(2,1));
//out(3,1)=(Rn*(1-e_squared)+x(3,1))*sin(x(1,1));

	double a = 6378137;
	double e_squared = .006694379990141316;

	double Rn = a/sqrt((1-e_squared*(sin(x(0,0))*sin(x(0,0)))));

	ublas::matrix<double> pos_ECEF(3,1);
	pos_ECEF(0,0) = (Rn + x(2,0))*cos(x(0,0))*cos(x(1,0));
	pos_ECEF(1,0) = (Rn + x(2,0))*cos(x(0,0))*sin(x(1,0));
	pos_ECEF(2,0) = (Rn*(1-e_squared)+x(2,0))*sin(x(0,0));

	return pos_ECEF;

}

ublas::matrix<double> GpsFunctions::attned2xyz(ublas::matrix<double> posECEF, ublas::matrix<double> attNED,double error)
{
	ublas::matrix<double> attECEF(3,1);
	attECEF.clear();

	double s1 = sin(attNED(0,0));
	double c1 = cos(attNED(0,0));
	double s2 = sin(attNED(1,0));
	double c2 = cos(attNED(1,0));
	double s3 = sin(attNED(2,0));
	double c3 = cos(attNED(2,0));

	ublas::matrix<double> CNEDBody(3,3);
	CNEDBody.clear();
	CNEDBody(0,0) = c2*c3;
	CNEDBody(0,1) = c2*s3;
	CNEDBody(0,2) = -1*s2;
	CNEDBody(1,0) = s1*s2*c3-c1*s3;
	CNEDBody(1,1) = s1*s2*s3+c1*c3;
	CNEDBody(1,2) = s1*c2;
	CNEDBody(2,0) = c1*s2*c3+s1*s3;
	CNEDBody(2,1) = c1*s2*s3-s1*c3;
	CNEDBody(2,2) = c1*c2;

	ublas::matrix<double> CECEFNED(3,3);
	CECEFNED.clear();
	CECEFNED = get_C_xyz2ned(posECEF,error);
	
	ublas::matrix<double> CECEFBODY(3,3);
	CECEFBODY.clear();
	CECEFBODY = prod(CNEDBody,CECEFNED);

	attECEF(0,0) = atan2(CECEFBODY(1,2),CECEFBODY(2,2));
	attECEF(1,0) = -asin(CECEFBODY(0,2));
	attECEF(2,0) = atan2(CECEFBODY(0,1),CECEFBODY(0,0));
	
	return attECEF;
}

ublas::matrix<double> GpsFunctions::attxyz2ned(ublas::matrix<double> posECEF, ublas::matrix<double> attECEF)
{

	/*
s1=sin(att_ecef(1,1));c1=cos(att_ecef(1,1));
s2=sin(att_ecef(2,1));c2=cos(att_ecef(2,1));
s3=sin(att_ecef(3,1));c3=cos(att_ecef(3,1));
C_ecef_body=[c2*c3 c2*s3 -s2; s1*s2*c3-c1*s3 s1*s2*s3+c1*c3 s1*c2; c1*s2*c3+s1*s3 c1*s2*s3-s1*c3 c1*c2];
C_ecef_ned=Get_C_ECEF_NED(pos_ecef,1);
C_ned_body=C_ecef_body*C_ecef_ned';

att_ned(1,1)=atan2(C_ned_body(2,3),C_ned_body(3,3));
att_ned(2,1)=-asin(C_ned_body(1,3));
att_ned(3,1)=atan2(C_ned_body(1,2),C_ned_body(1,1));*/
	double s1 = sin(attECEF(0,0));
	double s2 = sin(attECEF(1,0));
	double s3 = sin(attECEF(2,0));

	double c1 = cos(attECEF(0,0));
	double c2 = cos(attECEF(1,0));
	double c3 = cos(attECEF(2,0));

	ublas::matrix<double> C_ecef_body(3,3);
	C_ecef_body(0,0) = c2*c3;
	C_ecef_body(0,1) = c2*s3;
	C_ecef_body(0,2) = -s2;
	C_ecef_body(1,0) = s1*s2*c3 - c1*s3;
	C_ecef_body(1,1) = s1*s2*s3 + c1*c3;
	C_ecef_body(1,2) = s1*c2;
	C_ecef_body(2,0) = c1*s2*c3+s1*s3;
	C_ecef_body(2,1) = c1*s2*s3-s1*c3;
	C_ecef_body(2,2) = c1*c2;

	ublas::matrix<double> C_ecef_ned = get_C_xyz2ned(posECEF,.001);
	ublas::matrix<double> C_ned_body = prod(C_ecef_body,trans(C_ecef_ned));

	ublas::matrix<double> att_ned(3,1);
	att_ned(0,0) = atan2(C_ned_body(1,2),C_ned_body(2,2));
	att_ned(1,0) = -asin(C_ned_body(0,2));
	att_ned(2,0) = atan2(C_ned_body(0,1),C_ned_body(0,0));

	return att_ned;

}

double GpsFunctions::checkAngleWrap(double ang, double refAngleComp)
{
double outAng;
	outAng = ang;

	if(refAngleComp == 0)
	{
		if(ang > pi)
			outAng = ang-2*pi;

		if(ang < -pi)
			outAng = ang+2*pi;
	}

	if(refAngleComp == 1)
	{
		if(ang > pi/2)
			outAng = ang-pi;

		if(ang <(-pi/2))
			outAng = ang+pi;

	}
	
	return outAng;
}

double GpsFunctions::get_pseudorange_var(double c2n,int pseudorangeSel, int LSel)
{
	/*c2n=10^(c2n/10);


Fa=1;
Fb=1;
d=.5;
Bn=2;
lambda_c=293.05;

lambda_l1_L=.1903;
lambda_l2_L=.2442;
F=1;

%m
sigma_squared_atm=5.22^2;
%m/s
fe=3;

switch case_a

    case 1
        T=2;

        dll=lambda_c*sqrt(  ((4*Fa*(d^2)*Bn)/(c2n))  *  (2*(1-d)+((4*Fb*d)/(T*c2n)))  );

        var=dll^2+sigma_squared_atm;

    case 2
        T=5;
        switch case_b

            case 1

                fll=(lambda_l1_L/(2*pi*T))*sqrt( ((4*F*Bn)/c2n)*(1+(1/(T*c2n))));

                var=fll^2+((fe^2)/9);

            case 2

                fll=(lambda_l2_L/(2*pi*T))*sqrt( ((4*F*Bn)/c2n)*(1+(1/(T*c2n))));

                var=fll^2+((fe^2)/9);
        end
end*/

	c2n=pow(10,(c2n/10));


double Fa=1;
double Fb=1;
double d=.5;
double Bn=2;
double lambda_c=293.05;

double lambda_l1_L=.1903;
double lambda_l2_L=.2442;
double F=1;

//m
double sigma_squared_atm=pow(5.22,2);
//m/s
double fe=3;

double var = 1000;

if(pseudorangeSel == 0)
{

        double T=2;

        double dll=lambda_c*sqrt(  ((4*Fa*(pow(d,2))*Bn)/(c2n))  *  (2*(1-d)+((4*Fb*d)/(T*c2n)))  );

        var=pow(dll,2)+sigma_squared_atm;
}
else
if(pseudorangeSel == 1)
{
	double T=5;
       if(LSel == 0)
	   {

           //     fll=(lambda_l1_L/(2*pi*T))*sqrt( ((4*F*Bn)/c2n)*(1+(1/(T*c2n))));
	        double fll=(lambda_l1_L/(2*pi*T))*sqrt( ((4*F*Bn)/c2n)*(1+(1/(T*c2n))));

                var=pow(fll,2)+((pow(fe,2))/9);
	   }
	   else
	   if(LSel == 1)
	   {
                double fll=(lambda_l2_L/(2*pi*T))*sqrt( ((4*F*Bn)/c2n)*(1+(1/(T*c2n))));

                var=pow(fll,2)+((pow(fe,2))/9);
	   }
	   else
		   cout << "Bad parameter - get_pseudorange_var function, LSel parameter" << endl;
     
}
  else
		   cout << "Bad parameter - get_pseudorange_var function, pseudorangeSel parameter" << endl;

return var;

}

// int GpsFunctions::sv_pos_calc(MasEphemData Ephemerides, double svpos[], int sv, double gpsTime, double psr, double *satReturn) 
// {

//     /* Push ephemerides column 'sv' into properly named variables to
// facilitate calculation. */
//     double T_GD = Ephemerides.groupDelay[sv];
//     double t_oc = Ephemerides.clockCorrection[sv];
//     double a_f2 = Ephemerides.clockAging3[sv];
//     double a_f1 = Ephemerides.clockAging2[sv];
//     double a_f0 = Ephemerides.clockAging1[sv];
//     double C_rc = Ephemerides.cosOrbitRadius[sv];
//     double C_rs = Ephemerides.sinOrbitRadius[sv];
//     double C_uc = Ephemerides.cosLatitude[sv];
//     double C_us = Ephemerides.sinLatitude[sv];
//     double C_ic = Ephemerides.cosInclination[sv];
//     double C_is = Ephemerides.sinInclination[sv];
//     double Delta_n = Ephemerides.meanMotionDiff[sv];
//     double M_0 = Ephemerides.meanAnomaly[sv];
//     double e = Ephemerides.eccentricity[sv];
//     double sqrt_A = Ephemerides.semiMajorAxis[sv]; //MAKE SURE THIS PARAMETER IS INPUT AS THE SQRT!!!
//     double t_oe = Ephemerides.referenceTime[sv];
//     double Omega_0 = Ephemerides.rightAscension[sv];
//     double i_0 = Ephemerides.inclinationAngle[sv];
//     double omega = Ephemerides.perigeeArg[sv];
//     double dot_Omega = Ephemerides.ascensionRate[sv];
//     double Idot = Ephemerides.inclinationRate[sv];

//     //Check time of ephemeris, use latest one (?)
//     long toe = (long) t_oe; //double-check original intention
//     //toe can be replaced with NULL pointer check

//     //If ephemerides exist, toe will be a real number
//     if (toe == 0) {
//         cout << "GpsFunctions.cc - sv_pos_calc-> toe is NULL\r" << endl;
//         return (0);
//     }
//     //Calculate transit time (pseudorange/speed of light)
//     double transitTime = psr / C;
//     //Calculate transmit time (current time - transmit time)
//     double transmitTime = gpsTime - transitTime; // assumes time in seconds!!!
//     //----Find time difference--------------------------------------------
//     double dt = check_t(transmitTime - t_oc);
//     //--- Calculate clock correction--------------------------------------
//     double satClkCorr = (a_f2 * dt + a_f1) * dt + a_f0 - T_GD; //scope
//     double time = transmitTime - satClkCorr;
//     // Find satellite's position------------------------------------------
//     //Restore semi-major axis
//     double a = sqrt_A*sqrt_A;
//     //Time correction
//     double tk = check_t(time - t_oe);
//     //Initial mean motion
//     double n0 = sqrt(GM / pow(a, 3));
//     //Mean motion
//     double n = n0 + Delta_n;
//     //Mean anomaly
//     double M = M_0 + n * tk;
//     //Reduce mean anomaly to between 0 and 360 deg
//     M = fmod(M + 2 * gpsPi, 2 * gpsPi);
//     //Initial guess of eccentric anomaly
//     double E = M;
//     //----Iteratively compute eccentric anomaly---------------------------
//     int ii;
//     double E_old, dE, dtr;
//     for (ii = 0; ii < 10; ii++) {
//         E_old = E;
//         E = M + e * sin(E);
//         dE = fmod(E - E_old, 2 * gpsPi);
//         if (fabs(dE) < 1.e-12) {
//             //Necessary precision is reached, exit from the loop
//             break;
//         }
//     }
//     //Reduce eccentric anomaly to between 0 and 360 deg
//     E = fmod(E + 2 * gpsPi, 2 * gpsPi);
//     //Compute relativistic correction term
//     dtr = F*e; //*sqrt_A*sin(E);
//     //Calculate the true anomaly
//     double nu = atan2(sqrt(1 - e * e) * sin(E), cos(E) - e);
//     //Compute angle phi
//     double phi = nu + omega;
//     //Reduce phi to between 0 and 360 deg
//     phi = fmod(phi, 2 * gpsPi);
//     //Correct argument of latitude
//     double u = phi + C_uc * cos(2 * phi) + C_us * sin(2 * phi);
//     //Correct radius
//     double r = a * (1 - e * cos(E)) + C_rc * cos(2 * phi) + C_rs * sin(2 * phi);
//     //Correct inclination
//     double i = i_0 + Idot * tk + C_ic * cos(2 * phi) + C_is * sin(2 * phi);
//     //Compute the angle between the ascending node and the Greenwich meridian
//     double Omega = Omega_0 + (dot_Omega - Omegae_dot) * tk - Omegae_dot * t_oe - Omegae_dot*transitTime;
//     //Reduce to between 0 and 360 deg
//     Omega = fmod(Omega + 2 * gpsPi, 2 * gpsPi);
//     //----Compute satellite coordinates-----------------------------------
//     svpos[0] = cos(u) * r * cos(Omega) - sin(u) * r * cos(i) * sin(Omega); //ECEF X
//     svpos[1] = cos(u) * r * sin(Omega) + sin(u) * r * cos(i) * cos(Omega); //ECEF Y
//     svpos[2] = sin(u) * r * sin(i); //ECEF Z
//     //Include relativistic correction in clock correction-----------------
//     satClkCorr = (a_f2 * dt + a_f1) * dt + a_f0 - T_GD + dtr;
//     *satReturn = satClkCorr;

//     return 1;
// } //end sv_pos_calc

//calculates satellite position for Ephemerides boost matrix for one satellite
int GpsFunctions::sv_pos_calc(ublas::matrix<double> &Ephemerides, double svpos[], double svvel[], int sv, double gpsTime, double psr, double *satReturn) 
{
	//	[groupDelay, clockCorrection, clockAging3, clockAging2, clockAging1, cosOrbitRadius, sinOrbitRadius, cosLatitude, sinLatitude, cosInclination
	//	sinInclination, meanMotionDiff, meanAnomaly, eccentricity, semiMajorAxis, referenceTime, rightAscension, inclinationAngle, perigeeArg,
	//	ascensionRate, inclinationRate]'

    /* Push ephemerides column 'sv' into properly named variables to
facilitate gcalculation. */

    //cout << "Start sv_pos_calc.\n" << endl;

    double T_GD = Ephemerides(1,sv);
    double t_oc = Ephemerides(0,sv);
    double a_f2 = Ephemerides(20,sv);
    double a_f1 = Ephemerides(19,sv);
    double a_f0 = Ephemerides(18,sv);
    double C_rc = Ephemerides(10,sv);
    double C_rs = Ephemerides(11,sv);
    double C_uc = Ephemerides(8,sv);
    double C_us = Ephemerides(9,sv);
    double C_ic = Ephemerides(12,sv);
    double C_is = Ephemerides(13,sv);
    double Delta_n = Ephemerides(4,sv);
    double M_0 = Ephemerides(5,sv);
    double e = Ephemerides(6,sv);
    double sqrt_A = Ephemerides(3,sv); //MAsvE SURE THIS PARAMETER IS INPUT AS THE SQRT!!!
    double t_oe = Ephemerides(2,sv);
    double Omega_0 = Ephemerides(16,sv);
    double i_0 = Ephemerides(14,sv);
    double omega = Ephemerides(7,sv);
    double dot_Omega = Ephemerides(17,sv);
    double Idot = Ephemerides(15,sv);

    //cout << "Loaded ephemeris.\n" << endl;
    //Check time of ephemeris, use latest one (?)
    long toe = (long) t_oe; //double-check original intention
    //toe can be replaced with NULL pointer check

    //If ephemerides exist, toe will be a real number
    if (toe == 0) {
        cout << "GpsFunctions.cc - sv_pos_calc-> toe is NULL\r" << endl;
        return (0);
    }
    //Calculate transit time (pseudorange/speed of light)
    double transitTime = psr / C;
    //Calculate transmit time (current time - transmit time)
    double transmitTime = gpsTime - transitTime; // assumes time in seconds!!!
    //----Find time difference--------------------------------------------
    double dt = check_t(transmitTime - t_oc);
    //--- Calculate clock correction--------------------------------------
    double satClkCorr = (a_f2 * dt + a_f1) * dt + a_f0 - T_GD; //scope
    double time = transmitTime - satClkCorr;
    // Find satellite's position------------------------------------------
    //Restore semi-major axis
    double a = sqrt_A*sqrt_A;
    //Time correction
    double tk = check_t(time - t_oe);
    //Initial mean motion
    double n0 = sqrt(GM / pow(a, 3));
    //Mean motion
    double n = n0 + Delta_n;
    //Mean anomaly
    double M = M_0 + n * tk;
    //Reduce mean anomaly to between 0 and 360 deg
    M = fmod(M + 2 * gpsPi, 2 * gpsPi);
    //Initial guess of eccentric anomaly
    double E = M;
    //----Iteratively compute eccentric anomaly---------------------------
    int ii;
    double E_old, dE, dtr;
    for (ii = 0; ii < 10; ii++) {
        E_old = E;
        E = M + e * sin(E);
        dE = fmod(E - E_old, 2 * gpsPi);
        if (fabs(dE) < 1.e-12) {
            //Necessary precision is reached, exit from the loop
            break;
        }
    }
    //Reduce eccentric anomaly to between 0 and 360 deg
    E = fmod(E + 2 * gpsPi, 2 * gpsPi);
    //Compute relativistic correction term
    dtr = F*e; //*sqrt_A*sin(E);
    //Calculate the true anomaly
    double nu = atan2(sqrt(1 - e * e) * sin(E), cos(E) - e);
    //Compute angle phi
    double phi = nu + omega;
    //Reduce phi to between 0 and 360 deg
    phi = fmod(phi, 2 * gpsPi);
    //Correct argument of latitude
    double u = phi + C_uc * cos(2 * phi) + C_us * sin(2 * phi);
    //Correct radius
    double r = a * (1 - e * cos(E)) + C_rc * cos(2 * phi) + C_rs * sin(2 * phi);
    //Correct inclination
    double i = i_0 + Idot * tk + C_ic * cos(2 * phi) + C_is * sin(2 * phi);
    //Compute the angle between the ascending node and the Greenwich meridian
    double Omega = Omega_0 + (dot_Omega - Omegae_dot) * tk - Omegae_dot * t_oe - Omegae_dot*transitTime;
    //Reduce to between 0 and 360 deg
    Omega = fmod(Omega + 2 * gpsPi, 2 * gpsPi);

    double X = r*cos(u);
    double Y = r*sin(u);

    //----Compute satellite coordinates-----------------------------------
    svpos[0] = cos(u) * r * cos(Omega) - sin(u) * r * cos(i) * sin(Omega); //ECEF X
    svpos[1] = cos(u) * r * sin(Omega) + sin(u) * r * cos(i) * cos(Omega); //ECEF Y
    svpos[2] = sin(u) * r * sin(i); //ECEF Z
    //Include relativistic correction in clock correction-----------------
    satClkCorr = (a_f2 * dt + a_f1) * dt + a_f0 - T_GD + dtr;
    *satReturn = satClkCorr;

    double Edot = ( n0 + Delta_n )/( 1 - e*cos(E) );

    double phidot = ( sqrt(1- pow(e,2) )/( 1 - e*cos(E)) )*Edot;

    double udot = ( 1 + 2*C_us*cos(2*phi) - 2*C_uc*sin(2*phi) )*phidot;

    double rdot = 2*( C_rs*cos(2*phi) - C_rc*sin(2*phi) )*phidot + a*e*sin(E)*Edot;

    double idot = 2*( C_is*cos(2*phi) - C_ic*sin(2*phi) )*phidot + Idot;

    double Xdot = rdot*cos(u) - r*sin(u)*udot;

    double Ydot = rdot*sin(u) + r*cos(u)*udot;

    double Omegadot = dot_Omega - Omegae_dot;

    svvel[0] = Xdot*cos(Omega) - Ydot*cos(i)*sin(Omega) + Y*sin(i)*sin(Omega)*idot - svpos[1]*Omegadot;

    svvel[1] = Xdot*sin(Omega) + Ydot*cos(i)*cos(Omega) - Y*sin(i)*cos(Omega)*idot + svpos[0]*Omegadot;

    svvel[2] = Ydot*sin(i) + Y*cos(i)*idot;


    return 1;
} //end sv_pos_calc

UDmatrix GpsFunctions::GetUnitVectors(UDmatrix SatellitePositions, double EcefUserPosition[], UIvector Prns, int Obs) {
    UDmatrix UnitVectors(3, 33);
    double Range;
    double dx;
    double dy;
    double dz;

    UnitVectors.clear();
    //Is there a more efficient, compact way to do this?
    for (int k = 0; k < Obs; k++) {
        dx = SatellitePositions(0, Prns(k)) - EcefUserPosition[0];
        dy = SatellitePositions(1, Prns(k)) - EcefUserPosition[1];
        dz = SatellitePositions(2, Prns(k)) - EcefUserPosition[2];
        Range = pow(dx * dx + dy * dy + dz*dz, 0.5);

        UnitVectors(0, Prns(k)) = dx / Range;
        UnitVectors(1, Prns(k)) = dy / Range;
        UnitVectors(2, Prns(k)) = dz / Range;
    }
    return (UnitVectors);
}

UDvector GpsFunctions::LeastSquares(UDmatrix H, UDvector z) {

    int n=H.size2(); //Number of columns in H
    UDmatrix HtH_Inverse(n,n);
    UDvector x(n);
    InvertMatrix( (UDmatrix)prod(trans(H), H), HtH_Inverse);
    x=prod(HtH_Inverse, (UDvector)prod(trans(H), z));
//    cout << "z=" << z<<"\r"<< endl;
//    cout << "x=" << x<<"\r"<< endl;

    return(x);

}

UDvector GpsFunctions::WeightedLeastSquares(UDmatrix H, UDvector z, UDmatrix R)
{

//    xcarL1(:,kGps)=(H'*RcarL1_inv*H)^-1*H'*RcarL1_inv*zcarL1;

    int n=H.size2(); //Number of columns in H
    UDmatrix R_inv(R.size1(),R.size2());
    InvertMatrix(R, R_inv);
    UDmatrix HtH_Inverse(n,n);
    UDvector x(n);
    InvertMatrix( (UDmatrix)prod((UDmatrix)prod(trans(H), R_inv),H), HtH_Inverse);

//    cout << "Rinv:\r" << endl;
//    cout << R_inv << "\r" << endl << endl;
//
//    cout << "HtH_Inverse:\r" << endl;
//    cout << HtH_Inverse << "\r" << endl << endl;


    x=prod(HtH_Inverse, (UDvector)prod((UDmatrix)prod(trans(H), R_inv),z));
//    cout << "z=" << z<<"\r"<< endl;
//   cout << "x=" << x<<"\r"<< endl<< endl;

    return(x);

}

UIvector GpsFunctions::GetSignalList(UDvector RangeData, UDvector PhaseData, UDvector OldRangeData, UDvector OldPhaseData, double wavelength) {
    UIvector SignalsList(33);
    UIvector csigs(33);
    UIvector osigs(33);
    UIvector nsigs(33);
    UIvector noslipsigs(33);
    UDvector Nc(33);
    UDvector No(33);
    UDvector dN(33);

    Nc = RangeData - PhaseData;
    No = OldRangeData - OldPhaseData;
    dN = (Nc - No) / wavelength;

//    cout << "dN:\r"<< endl;
//    cout << dN << "\r"<< endl;

    double tol = 1.0;
    double negtol = -1.0;
    csigs(0) = 0;
    osigs(0) = 0;
    noslipsigs(0) = 0;
//    for (int k = 1; k < 33; k++) {
//        //Find current list of signals
//        if (!isnan(RangeData(k))) {
//            csigs(k) = 1;
//        } else {
//            csigs(k) = 0;
//        }
//        //Find old list of signals
//        if (isnan(OldRangeData(k))) {
//            osigs(k) = 1;
//        } else {
//            osigs(k) = 0;
//        }
//        //Check for cycle slips
//        if (dN(k)>-tol && dN(k) < tol) {
//            noslipsigs(k) = 1;
//        } else {
//            noslipsigs(k) = 0;
//        }
//    }
//    nsigs = element_prod(csigs, osigs);
//    SignalsList = nsigs + noslipsigs;
    for (int k = 1; k < 33; k++) {
        //Find current list of signals
        if (!boost::math::isnan(RangeData(k))) {
            csigs(k) = 1;
        } else {
            csigs(k) = 0;
        }
        //Find old list of signals
        if (!boost::math::isnan(OldRangeData(k))) {
            osigs(k) = 1;
        } else {
            osigs(k) = 0;
        }
        //Check for cycle slips
        if (dN(k)>-tol && dN(k) < tol) {
            noslipsigs(k) = 1;
        } else {
            noslipsigs(k) = 0;
        }
    }

    nsigs = element_prod(csigs, osigs);
    SignalsList = element_prod(nsigs, noslipsigs);
//    cout <<"csigs, osigs, nsigs\r"<< endl;
//    cout << csigs << "\r" << endl;
//    cout << osigs << "\r" << endl;
//    cout << nsigs << "\r" << endl;
//    cout << noslipsigs << "\r" << endl;

    return (SignalsList);

}

// UDmatrix GpsFunctions::GetSatellitePosition(MasEphemData ephemerides, UDvector Pseudorange, UIvector Prns, int Obs, double GpsTime) {
//     UDmatrix SatellitePositions(3, 33);
//     double svpos[3] = {0.};
//     double svClockCorrection;

//     SatellitePositions.clear();
//     for (int k = 0; k < Obs; k++) {
//         sv_pos_calc(ephemerides, svpos, Prns(k), GpsTime, Pseudorange(Prns(k)), &svClockCorrection);
//         //Put satellite positions into matrix
//         for (int j = 0; j < 3; j++) {
//             SatellitePositions(j, Prns(k)) = svpos[j];
//         }
//         //	cout << "Prns= " << Prns(k) << "\tSatPos(1)= " << SatellitePositions(0,Prns(k)) << "\tSatPos(2)= " <<SatellitePositions(1,Prns(k)) <<  "\tSatPos(1)= " << SatellitePositions(2,Prns(k)) << "\r" << endl;
//     }
//     return (SatellitePositions);
// }

UDmatrix GpsFunctions::GetSatellitePosition(ublas::matrix<double> ephem_mat, std::vector<double> &Pseudorange, UIvector Prns, int Obs, double GpsTime) {
    
    //cout << "Start satellite position calculation.\n" << endl;
    UDmatrix SatellitePositions(7, 32);
    double svpos[3] = {0.};
    double svvel[3] = {0.};
    double svClockCorrection;

    //cout << Pseudorange[0] << endl;
    //cout << Pseudorange[1] << endl;

    SatellitePositions.clear();
    for (int k = 0; k < Obs; k++) {

        //cout << "Start calculation.\n" << endl;

        //cout << "prns = " << Prns(k) - 1 << endl;
        //cout << "pseudorange = " << Pseudorange[Prns(k)-1] << endl;
        sv_pos_calc( ephem_mat, svpos, svvel, Prns(k) - 1, GpsTime, Pseudorange[Prns(k)-1], &svClockCorrection);
        //cout << "Position calculated.\n" << endl;
        //Put satellite positions into matrix
        for (int j = 0; j < 3; j++) {
            SatellitePositions(j, Prns(k) - 1) = svpos[j];
            SatellitePositions(j+4, Prns(k) - 1) = svvel[j];
        }
        SatellitePositions(3,Prns(k) - 1) = svClockCorrection;
        //cout << "Prns= " << Prns(k) << "\tSatPos(1)= " << SatellitePositions(0,Prns(k)-1) << "\tSatPos(2)= " <<SatellitePositions(1,Prns(k)-1) <<  "\tSatPos(3)= " << SatellitePositions(2,Prns(k)-1) << "\r" << endl;
    }
    return (SatellitePositions);
}


UDvector GpsFunctions::ExpandRangeData(double RangeData[14], int FreshDataPrns[14]) {
    UDvector ExpandedData(33);
    //Initialize all elements to NaN
    for (int k = 0; k < 33; k++) {
        ExpandedData(k) = std::numeric_limits<double>::quiet_NaN();
    }
    //Extract range data and place in row/colum associated with the PRN of the data
    for (int k = 0; k < 14; k++) {
        for (int j = 1; j < 33; j++) {
            if (FreshDataPrns[k] == j) {
                ExpandedData(j) = RangeData[k];
                break;
            }
        }
    }
    return (ExpandedData);
}
//void GenerateEphemerisDatabase(MasEphemData Ephemerides);
// void GpsFunctions::GenerateEphemerisDatabase(MasEphemData EphemerisDatabase, MasEphemData Ephemerides){
//     	for(int sv=1;sv<33;sv++){
// 		//Need to add a check on the GPS week number, or else new ephemerides
// 		//won't get loaded in the event of a week rollover
    
// 		if(Ephemerides.referenceTime[sv]>EphemerisDatabase.referenceTime[sv]){
// 			EphemerisDatabase.groupDelay[sv]=Ephemerides.groupDelay[sv];
//             EphemerisDatabase.clockCorrection[sv]=Ephemerides.clockCorrection[sv];
//             EphemerisDatabase.clockAging3[sv]=Ephemerides.clockAging3[sv];
//             EphemerisDatabase.clockAging2[sv]=Ephemerides.clockAging2[sv];
// 			EphemerisDatabase.clockAging1[sv]=Ephemerides.clockAging1[sv];
// 			EphemerisDatabase.cosOrbitRadius[sv]=Ephemerides.cosOrbitRadius[sv];
// 			EphemerisDatabase.sinOrbitRadius[sv]=Ephemerides.sinOrbitRadius[sv];
// 			EphemerisDatabase.cosLatitude[sv]=Ephemerides.cosLatitude[sv];
// 			EphemerisDatabase.sinLatitude[sv]=Ephemerides.sinLatitude[sv];
// 			EphemerisDatabase.cosInclination[sv]=Ephemerides.cosInclination[sv];
// 			EphemerisDatabase.sinInclination[sv]=Ephemerides.sinInclination[sv];
// 			EphemerisDatabase.meanMotionDiff[sv]=Ephemerides.meanMotionDiff[sv];
// 			EphemerisDatabase.meanAnomaly[sv]=Ephemerides.meanAnomaly[sv];
// 			EphemerisDatabase.eccentricity[sv]=Ephemerides.eccentricity[sv];
// 			EphemerisDatabase.semiMajorAxis[sv]=Ephemerides.semiMajorAxis[sv];
// 			EphemerisDatabase.referenceTime[sv]=Ephemerides.referenceTime[sv];
// 			EphemerisDatabase.rightAscension[sv]=Ephemerides.rightAscension[sv];
// 			EphemerisDatabase.inclinationAngle[sv]=Ephemerides.inclinationAngle[sv];
// 			EphemerisDatabase.perigeeArg[sv]=Ephemerides.perigeeArg[sv];
// 			EphemerisDatabase.ascensionRate[sv]=Ephemerides.ascensionRate[sv];
// 			EphemerisDatabase.inclinationRate[sv]=Ephemerides.inclinationRate[sv];
// 			//EphemerisDatabase.timeOfWeek[sv]=Ephemerides.timeOfWeek[sv];
// 		}
// 	}
        
// }

UDmatrix GpsFunctions::GetGeometryMatrix(UDmatrix UnitVectors, UIvector Prns, int Obs) {
    UDmatrix GeometryMatrix(Obs, 3);

    GeometryMatrix.clear();
    for (int k = 0; k < Obs; k++) {
        GeometryMatrix(k, 0) = UnitVectors(0, Prns(k));
        GeometryMatrix(k, 1) = UnitVectors(1, Prns(k));
        GeometryMatrix(k, 2) = UnitVectors(2, Prns(k));
        //GeometryMatrix(k, 3) = 1;
    }
    return (GeometryMatrix);
}

UDvector GpsFunctions::GetMeasurementVector(UDvector carL1, UDvector old_carL1, UDmatrix svpos, UDmatrix oldsvpos, UDmatrix geometry_matrix, UIvector Prns, int Obs) {

    UDvector z(Obs);
    double dcar;
    double dsvpos;
//    cout << "obs=" << Obs << endl;
//    cout << Prns << "\r" << endl;
//    cout << carL1 << "\r" << endl;
//    cout << old_carL1 << "\r" << endl << endl;
//
//    cout << "svpos:" << "\r" << endl;
//    cout << svpos << "\r" << endl;
//    cout << "oldsvpos:" << "\r" << endl;
//    cout << oldsvpos << "\r" << endl << endl;
//
//    cout << geometry_matrix << "\r" << endl << endl;
    for (int k = 0; k < Obs; k++) {

//        cout << carL1(Prns(k)) << "\t" << old_carL1(Prns(k)) << "\r" << endl;
        dcar = carL1(Prns(k)) - old_carL1(Prns(k));  //Time difference of carrier
        dsvpos = 0.; //Holds time difference of sv position
        for (int j = 0; j < 3; j++) {
//            cout << j << "\t"<< k << "\t";
//            cout <<geometry_matrix(k,j) << "\t";
//            cout << svpos(j, Prns(k)) << "\t";
//            cout <<oldsvpos(j, Prns(k)) <<"\t";
            dsvpos = dsvpos + geometry_matrix(k, j)*(svpos(j, Prns(k)) - oldsvpos(j, Prns(k)));
//            cout << dsvpos << "\r" << endl;
        }
        z(k) = dcar - dsvpos;
//        cout << k << "\tz(k)=" << z(k) << "\tdcar=" << dcar<< "\tdsvpos=" << dsvpos <<"\r"<< endl;
    }
    return (z);
}

UDmatrix GpsFunctions::GetMeasurementCovariance(int Obs, UIvector Prns, UDvector Rcvr1CnoL1) {

    double wo = 22.94455;
    double Bn = 0.7845 * wo; //%Third order loop filter  (18Hz)
    double T = 1 / 200.; //%Predetection integration time
    int F1 = 1; //%Discriminator correlator factor
    int F2 = 1; //%Discriminator type factor
    double d = 0.5; //%Correlator spaceing (chips)
    int Bnc = 2; //%Code loop noise bandwidth (Hz)
    double ccr = 293.05; //%Code chipping rate (chips/s)
    double var_atm=0.;
    double var_car_atm=0.;

    double r1CNoL1;
    double waveL1=0.190293672798365;
    double waveL2=0.244210213424568;

    UDmatrix R(Obs, Obs);
    R.clear();

    for (int n = 0; n < Obs; n++) {
        r1CNoL1 = pow(10, (Rcvr1CnoL1(Prns(n)) / 10));

        //L1 carrier variance (m^2)
        R(n,n) = var_car_atm + pow((waveL1 / (2 * gpsPi)),2)*(Bn / r1CNoL1)*(1 + (1 / (2 * T * r1CNoL1)));
    }

    return (R);
}

// UIvector GpsFunctions::GetNumEphems(MasEphemData ephemerides) {
//     UIvector EphemList(33);
//     EphemList(0) = 0;
//     for (int k = 1; k < 33; k++) {
//         if (ephemerides.semiMajorAxis[k]) {
//             EphemList(k) = 1;
//         } else {
//             EphemList(k) = 0;
//         }
//     }
//     return (EphemList);
// }

UIvector GpsFunctions::GetPrns(UIvector sigs, int obs) {
    UIvector prns(obs);
    int j = 0;
    for (int k = 0; k < 33; k++) {
        if (sigs(k)) {
            prns(j) = k;
            j++;
        }
    }
    return (prns);
}

UDmatrix GpsFunctions::skew(UDvector in_vec) {

    UDmatrix out_mat(3,3);
    out_mat.clear();

    // put in_vec into skew symmetric matrix form
    out_mat(0,0) = 0.;
    out_mat(0,1) = -in_vec(2);
    out_mat(0,2) = in_vec(1);
    out_mat(1,0) = in_vec(2);
    out_mat(1,1) = 0.;
    out_mat(1,2) = -in_vec(0);
    out_mat(2,0) = -in_vec(1);
    out_mat(2,1) = in_vec(0);
    out_mat(2,2) = 0.;

    return (out_mat);
}





//
// File:   GpsFunctions.h
// Author: gavlab
//
// Created on March 18, 2009, 1:35 PM
//

#ifndef _GPSFUNCTIONS_H
#define	_GPSFUNCTIONS_H

namespace GpsFunctions {
}; //tells it its going to use it
//struct MasEphemData;

#include <iostream>
#include <cmath>
#include <float.h>
#include <fstream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/function.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/assert.hpp>
//#include "gpsstructs.h"

#include <boost/assign/std/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "InvertMatrix.h"


using namespace boost::numeric::ublas;
using namespace std;
using namespace GpsFunctions;
using namespace boost::posix_time;
using namespace boost::gregorian;

typedef boost::numeric::ublas::vector<double> UDvector;
typedef boost::numeric::ublas::matrix<double> UDmatrix;
typedef boost::numeric::ublas::vector<int> UIvector;
typedef boost::numeric::ublas::matrix<int> UImatrix;
typedef boost::numeric::ublas::identity_matrix<double> IDmatrix;
typedef boost::numeric::ublas::triangular_matrix<double, boost::numeric::ublas::upper, boost::numeric::ublas::row_major> UTriMatrix;


namespace ublas = boost::numeric::ublas;


static const double eRad = 6378137; //WGS-84 Equatorial radius (m)
static const double ecc = .08181919;//eccentricity 
static const double pi = boost::math::constants::pi<double>();//boost pi

/*template<class T>
 int InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse)
 {
	 typedef permutation_matrix<std::size_t> pmatrix;
 	//using namespace boost::numeric::ublas;
 	
        //check to make sure the matrix is square and invertible
        if (input.size1() != input.size2())
            return -21;
               
 	// create a working copy of the input
 	matrix<T> A(input);
        // matrix<T> A;
        // A.assign(input);
 	// create a permutation matrix for the LU-factorization
 	pmatrix pm(A.size1());

        // std::cout << "invertmatrix" << std::endl;
        //std::cout << A << std::endl;
 	// perform LU-factorization
	int res = lu_factorize(A, pm);	
 	if( res != 0 ) return -22;

 	// create identity matrix of "inverse"
 	inverse.assign(ublas::identity_matrix<T>(A.size1()));

 	// backsubstitute to get the inverse
 	lu_substitute(A, pm, inverse);

        return 0;
 }*/

namespace GpsFunctions {
#ifndef PI
#define PI 3.14159265358979 // pi
#endif

    // used to convert lat and long to UTM coordinates
#define GRAD_A_RAD(g) ((g)*0.0174532925199433)

    double checkWrap(double angle);
    int timeStamp(int &gps_weeks, int &gps_ms);

    const double gpsPi = 3.1415926535898; // Pi used in the GPS coordinate system
    const double Omegae_dot = 7.2921151467e-5; // Earth rotation rate, [rad/s]
    const double GM = 3.986005e14; // Earth's universal [m^3/s^2]
    const double F = -4.442807633e-10; // Constant, [sec/(meter)^(1/2)]
    const double C = 299792458; //Speed of light... that's fast!




    // Some constants used by these functions.
    static const double fe = 500000.0;
    static const double ok = 0.9996;

    const double FOURTHPI(PI / 4);
    const double deg2rad(PI / 180.0);
    const double rad2deg(180.0 / PI);


    // An array containing each vertical UTM zone.
   // static char cArray[] = "CDEFGHJKLMNPQRSTUVWX";

    double check_t(double time);

    double wrap_heading(double heading);

    int wgslla2utm(double Lat, double Long, double &northing, double &easting, int &zone, bool &north);
    int wgslla2xyz(double xyz[], double wlat, double wlon, double walt);
    int rot(double R[3][3], double angle, int axis);
    int rot(ublas::matrix<double> &R, double angle, int axis);
    UDmatrix skew(UDvector in_vec);
    int wgsxyz2enu(double enu[], double xyz[], double reflat, double reflon, double refalt);
    int wgsdiffxyz2diffenu(double diffenu[], double diffxyz[], double reflat, double reflon);
    int wgsdiffxyz2diffenu(ublas::vector<double> &diffenu, std::vector<double> &diffxyz, double reflat, double reflon);
    int wgslla2enu(double enu[], double lat, double lon, double alt, double reflat, double reflon, double refalt);
//    int wgsxyz2lla(double &lat, double &lon, double &alt, double xyz[]);
    int wgsxyz2lla(double &lat, double &lon, double &alt, std::vector<double> &xyz);

	//boost implemented functions (NED rather than ENU)
	ublas::matrix<double> get_C_xyz2ned(ublas::matrix<double> x,double error);
	ublas::matrix<double> get_C_wgslla2ned(ublas::matrix<double> x);
	ublas::matrix<double> xyz2wgslla(ublas::matrix<double> x,double error);
	ublas::matrix<double> wgslla2xyz(ublas::matrix<double> x);
	ublas::matrix<double> attned2xyz(ublas::matrix<double> posECEF, ublas::matrix<double> attNED,double error);
	ublas::matrix<double> attxyz2ned(ublas::matrix<double> posECEF, ublas::matrix<double> attECEF);	

	///////////////////////////
	//special functions
	//////////////////////////
	//if option=0 +/- pi
//%if option=1 +/- pi/2
	double checkAngleWrap(double ang, double refAngleComp);

	//get pseudorange variance and pseuodrange_dot variance
//case_a 0:pseudorange
//       1:pseudorange_dot
//if pseudorange_dot
//       0:L1 frequency
//       1:L2 frequency
	double get_pseudorange_var(double c2n,int pseudorangeSel, int LSel);

	//calculate sv_pos_calc
	//returns 0 if ephemerides don't exist
	//1 if all went smoothly
//	int sv_pos_calc(MasEphemData Ephemerides, double svpos[], int sv, double gpsTime, double psr, double *satReturn);
	int sv_pos_calc(ublas::matrix<double> &Ephemerides, double svpos[], double svvel[], int sv, double gpsTime, double psr, double *satReturn);
	UDmatrix GetUnitVectors(UDmatrix SatellitePositions, double EcefUserPosition[], UIvector Prns, int Obs);
	UDvector LeastSquares(UDmatrix H, UDvector z);
	UDvector WeightedLeastSquares(UDmatrix H, UDvector z, UDmatrix R);
	UIvector GetSignalList(UDvector RangeData, UDvector PhaseData, UDvector OldRangeData, UDvector OldPhaseData, double wavelength);
//	UDmatrix GetSatellitePosition(MasEphemData Ephemerides, UDvector Pseudorange, UIvector Prns, int Obs, double GpsTime);
	UDmatrix GetSatellitePosition(ublas::matrix<double> ephem_mat, std::vector<double> &Pseudorange, UIvector Prns, int Obs, double GpsTime);

	UDvector ExpandRangeData(double RangeData[14], int FreshDataPrns[14]);
//    void GenerateEphemerisDatabase(MasEphemData EphemerisDatabase, MasEphemData Ephemerides);
    UDmatrix GetGeometryMatrix(UDmatrix UnitVectors, UIvector Prns, int Obs);
    UDvector GetMeasurementVector(UDvector carL1, UDvector old_carL1, UDmatrix svpos, UDmatrix oldsvpos, UDmatrix geometry_matrix,UIvector Prns, int Obs);
    UDmatrix GetMeasurementCovariance(int Obs, UIvector Prns, UDvector Rcvr1CnoL1);
//    UIvector GetNumEphems(MasEphemData ephemerides);
    UIvector GetPrns(UIvector sigs, int obs);
}
#endif	/* _GPSFUNCTIONS_H */


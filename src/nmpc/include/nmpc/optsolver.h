// OptSolver.h

//Some Visual Studio dll stuff
#ifndef __OPTSOLVER_H_INCLUDED__   // if x.h hasn't been included yet...
#define __OPTSOLVER_H_INCLUDED__   //   #define this so the compiler knows it has been included



//#include <QtCore/QtGlobal>

// #if defined(OPTSOLVER_LIBRARY)
// #  define OPTSOLVER_EXPORT Q_DECL_EXPORT
// #else
// #  define OPTSOLVER_EXPORT Q_DECL_IMPORT
// #endif


//More Visual Studio dll stuff
//#ifdef OptSolver_EXPORTS
//#define optFunc __declspec(dllexport)
//#else
//#define optFunc __declspec(dllimport)
//#endif

#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <eigen3/Eigen/Dense> //seems to be fixed...
//#include "gOptPathUGV.h"
using Eigen::MatrixXd;

using namespace std;
using namespace Eigen;


namespace WarEagle
{
    class OptSolver
    {
    public:
        //Vehicle Parameters
        double m; //kg
        double hcg; //m
        double Ixx; //kg*m^2
        double cfRoll;
        double k_phiAux; // Nm/rad
        double a; //m
        double b; //m
        double iPsi; //kg*m^2

        //Magic Formula Parameters
        double B;
        double C;
        double D;


        int itNum;

        //Target and obstacle locations
        //double target(2);
        //double obstacle(2);

        //Control Horizon
        //double tLook;

        //Initial guesses for control input parameters
        //vector <double> u(4);
        vector <double> patSearchEick(vector <double>& x0, vector <double>& x0_begin, vector <double>& u, double tLook, vector <double>& target, vector <double>& obstacle, double uLast, double Vmax, double Vmin);

        OptSolver();
        ~OptSolver();

        // These are all 7 component vectors
        vector <double> x0Dot;
        vector <double> x0_k2Vec;
        vector <double> x0_k3Vec;
        vector <double> x0_k4Vec;  
        vector <double> k1;
        vector <double> k2;
        vector <double> k3;
        vector <double> k4;
        vector <double> xNext_cc; //vector in cost calc function
        vector <double> xNextVec;

        
        vector <double> costVec; // 2 components

        vector <double> uu; // 4 components
        vector <double> uNew; // 4 components
        vector <double> costVecNew; // 2 components
        vector <double> costVecPat; // 2 components
        vector <double> uOpt; // 4 components

        //Changing some stuff up for obstacle detection
        vector <double> xObVec;
        vector <double> yObVec;
        vector <double> distObVec;

        //And now the EigenMatrixLib stuff
        VectorXd xVecAug;
        MatrixXd Q;
        MatrixXd pat;
        MatrixXd uPat;

        double dt;








    private:
        double Vx;
        double tt;
        double deltaF;
        double dF1;
        double dF2;
        double uIn[2];
        double dF;
        double x;
        double y;
        double phi;
        double phiDot;
        double psi;
        double psiDot;
        double Vy;
        double alphaF;
        double alphaR;
        double Fyf;
        double Fyr;
        double psiDotDot;
        double VyDot;
        double betaDot;
        double yDot;
        double xDot;
        double phiDotDot;
        double x0DotArr[7];

        int k;

        double x1;
        double x2;
        double x3;
        double x4;
        double x5;
        double x6;
        double x7;
        double x0_k2[7];
        double x0_k3[7];
        double x0_k4[7];

        double kGroup[7];
        double xNext[7];

        // These guys can get assigned in the constructor
        double rough;
        double grade;
        double Vmax;
        double Vmin;
        double alpha;
        double beta;
        double gamma;
        double eta;
        double epsilon;

        int rwCnt;
        int colCnt;

        double t;
        double xprime;
        double yprime;
        double xOb;
        double yOb;
        double distOb;
        double distTarg;
        double J;
        double Jtot;
        int conVio;

        double delThresh;
        double del_uThresh;
        double delta_uThresh;
        double rollThresh;
        double rAvoid; 

        double theta;
        int lambda;
        double del_u;
        double costLow;
        int success;

        int p;
        int n;

        bool obDetected; //boolean for if an obstacle has been detected in field of view
        double minDistOb; //minimum distance to a detected obstacle




        vector <double> nonlinVehMod_rk4(double t, vector <double>& x0, int k, vector <double>& u, double tLook);

        vector <double> rk4(vector <double>& x0, double t, vector <double>& u, double tLook);

        vector <double> costCalc(vector <double>& x0, vector <double>& x0_begin, vector <double>& u, double tLook, double rollThresh, double rAvoid, vector <double>& target, vector <double>& obstacle, double uLast, double Vmax, double Vmin);

        // vector <double> patSearchEick(double dt, vector <double>& x0, vector <double>& x0_begin, vector <double>& u, double tLook, vector <double>& target, vector <double>& obstacle, double uLast);






    };
}


#endif // __MYCLASS_H_INCLUDED__

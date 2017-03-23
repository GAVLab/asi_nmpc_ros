#ifndef EKF_H
#define EKF_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/function.hpp>
#include <iostream>
#include <fstream>
#include "matrix2.h"

namespace ublas = boost::numeric::ublas;
using namespace std;

//////////////////////////////
//EKF class for extended Kalman filter
// - assumes nonlinear time update approximated with 4th order runge kutta
// - assumes linear measurement update function (have to do more for nonlinear version)
///////////////////////////////
class EKF
{
public:
	//Functions
	EKF();
	EKF(ublas::matrix <double> initX,ublas::matrix <double> initP,ublas::matrix<double> initR,ublas::matrix<double> initQ);
	~EKF();
	void reset();
	void init(ublas::matrix <double> initX,ublas::matrix <double> initP);
	//void initTuningStates()
	//{
		///////////////////////
		//add code for selecting tuning states / matrix sizes here
		//////////////////////
		//Q.resize(numStates,numStates,false);
		//Q.clear();

		//R.resize(numStates,numStates,false);
		//R.clear();

		//P.resize(numStates,numStates,false);
		//P.clear();

		//R(0,0) = .001;

		//Q(0,0) = .001;

		//P(0,0) = .01;
	//}
	//void initStates(ublas::matrix <double> initX);

	////////////////////////////
	//measEqu
	//mechanization equation for nonlinear time update
	////////////////////////////
	//ublas::matrix<double> mechEqu(ublas::matrix <double> u,ublas::matrix <double> x)
	//{
	//	//ublas::matrix<double> returnVar;
	//	//returnVar.resize(numStates,1,false);

	//	////////////////////////////////////////////////
	//	////Edit Mechanization equations here
	//	////////////////////////////////////////////////
	//	//returnVar(0,0) = u(0,0);
	//	////returnVar(1,0) = u(2,0) - xhat(2,0);
	//	////returnVar(2,0) = 0;
	//	////returnVar(3,0) = u(1,0) - xhat(4,0);
	//	////returnVar(4,0) = 0;

	//	return returnVar;

	//}

	////////////////////////////
	//measMechEqu
	//mechanization equation for nonlinear measurement update
	////////////////////////////
	//ublas::matrix<double> measMechEqu(ublas::matrix <double> x)
	//{
		//ublas::matrix<double> returnVar;
		//returnVar.resize(numStates,1,false);

		////////////////////////////////////////////////
		////Edit Mechanization equations here
		////////////////////////////////////////////////
		//returnVar(0,0) = 0;
		////returnVar(1,0) = u(2,0) - xhat(2,0);
		////returnVar(2,0) = 0;
		////returnVar(3,0) = u(1,0) - xhat(4,0);
		////returnVar(4,0) = 0;

		//return returnVar;

	//}


	//ublas::matrix<double> getA(ublas::matrix<double> u, ublas::matrix<double> xhat, double dt)
	//{
	//	ublas::matrix<double> returnVar;
	//	returnVar.resize(numStates,numStates,false);
	//	returnVar.clear();

	//	////////////////////////////////////
	//	//Edit A matrix (Jacobian of equations of motion) or otherwise
	//	////////////////////////////////////
	//	returnVar(0,0) = 1;
	//	//returnVar(1,1) = 1;
	//	//returnVar(2,2) = 1;
	//	//returnVar(3,3) = 1;
	//	//returnVar(4,4) = 1;		

	//	//returnVar(0,1) = dt*sin(xhat(3,0));
	//	//returnVar(0,3) = xhat(1,0)*cos(xhat(3,0));
	//	//returnVar(1,2) = -dt;
	//	//returnVar(3,4) = -dt;
	//	////////////////////////////////////

	//	return returnVar;
	//}

	void timeUpdate(boost::function3<ublas::matrix<double>,ublas::matrix<double>, ublas::matrix<double>, double>, ublas::matrix<double> u, ublas::matrix<double> A, ublas::matrix<double> W, ublas::matrix<double> Q, double dt);
	void measUpdate(boost::function1<ublas::matrix<double>,ublas::matrix<double> >,ublas::matrix<double> z,ublas::matrix<double> H, ublas::matrix<double> V, ublas::matrix<double> R);
	//void timeUpdate(ublas::matrix<double>(*)(ublas::matrix<double> x,ublas::matrix<double> u,double dt), ublas::matrix<double> u, ublas::matrix<double> A, ublas::matrix<double> W, ublas::matrix<double> Q, double dt);
	//void measUpdate(ublas::matrix<double>(*)(ublas::matrix<double> x),ublas::matrix<double> z,ublas::matrix<double> H, ublas::matrix<double> V, ublas::matrix<double> R);
	
	
	void writeColumnMatrixToFile(const char* filename,const char* matrixName,ublas::matrix<double> matrix);
	void writeMatrixToFile(const char* filename,const char* matrixName,ublas::matrix<double> matrix);
	ublas::matrix<double> loadData(const char* filename);
	void setOptions(bool recordAll, bool eraseAll);
	void writeDoubleToFile(const char* filename,double value);
	void setP(ublas::matrix<double> Pp);
	void setXhat(ublas::matrix<double> xhatp);
	ublas::matrix<double> getP();
	ublas::matrix<double> getXhat();
	//void setA(ublas::matrix<double> AParam);

	//boost::function3<ublas::matrix<double>,ublas::matrix<double>, ublas::matrix<double>, double> f;
	//boost::function1<ublas::matrix<double>,ublas::matrix<double>> h;


	//Variables

	//time update
//	ublas::matrix<double> A; // n x n state transition matrix / Jacobian matrix of partial derivatives of f with respect to x
//	ublas::matrix<double> W;// Jacobian matrix of partial derivatives of h with respect to x
//	ublas::matrix<double> B; // n x l control input transition matrix
//	ublas::matrix<double> H; // m x n measurement transition matrix / Jacobian matrix of partial derivatives of h with respect to x
//	ublas::matrix<double> u;// current inputs
//	ublas::matrix<double> uPrev;//previous inputs (for rungekutta) at time k-1
//	ublas::matrix<double> Q;//process noise covariance
//	ublas::matrix<double> R;//measurement noise covariance
	//	 - R is n x n with measurement covariances corresponding to each state
//	ublas::matrix<double> K;// n x m gain / blending factor

//	ublas::matrix<double> V;// Jacobian matrix of parital derivatives of h with respect to v

	int n;//number of states
	int m;//number of measurements
	int l;// number of inputs

private:
	//Functions
	ublas::matrix<double> inverse(ublas::matrix<double> matrix);
	//ublas::matrix<double> rungeKuttaTime(ublas::matrix<double> u,ublas::matrix<double> x,ublas::matrix<double>(*)(ublas::matrix<double> u, ublas::matrix<double> x),double dt);
	//ublas::matrix<double> rungeKuttaMeas(ublas::matrix<double> x,ublas::matrix<double>(*h)(ublas::matrix<double> x));
	//ublas::matrix<double> getH(int numMeas,ublas::matrix <bool> measStateSel);
	
	//ublas::matrix<double> getH(ublas::matrix<double> meas)
	//{
	//	int numMeas = meas.size1();
	//	ublas::matrix<double> returnVar;
	//	returnVar.resize(numMeas,numStates,false);
	//	returnVar.clear();

	//	////////////////////////////////////
	//	//Edit H matrix (measurement nonlinear function)
	//	////////////////////////////////////
	//	returnVar(0,0) = 1;
	//	////////////////////////////////////

	//	return returnVar;	
	//}

	//////////////////////////////////
	//recordAll - records variables of the filter in a text file
	//linearTime - use the linear time update
	//linearMeas - use the linear measurement update


	void recordPXhat();// P and xhat
	void recordTimeData(ublas::matrix<double> u, ublas::matrix<double> A, ublas::matrix<double> W, ublas::matrix<double> Q, double dt);//variables exclusive to time update
	void recordMeasurementData(ublas::matrix<double> z,ublas::matrix<double> H,ublas::matrix<double> V,ublas::matrix<double> R);//variables exclusive to measurement update
	void eraseData();
	void clearFile(const char* filename);


	//Private Variables
	ublas::matrix<double> xhat; //state vector n x 1
	ublas::matrix<double> P;//error covariance
	int counter; // iteration counter - increments on time and measurement update


	bool record;
	//bool timeLinearity;
	//bool measLinearity;
};

#endif
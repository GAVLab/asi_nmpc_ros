#include "EKF.h"

/////////////////////////////////
//Constructor
//initialize with 1 state
/////////////////////////////////
EKF::EKF()
{
	cout << "EKF Constructor" << endl;

	n = 1;
	m = 1;

	xhat.resize(n,1,false);
	xhat.clear();

	P.resize(n,n,false);
	P.clear();

	//uPrev.resize(n,1,false);
	//uPrev.clear();

	//initTuningStates();

	counter = 0;

	record = true;

}

/////////////////////////////////
/*
Normal Constructor
initX
- initX is n x 1 column matrix
- each index is filled with initial value
*/
/////////////////////////////////
EKF::EKF(ublas::matrix <double> initX, ublas::matrix <double> initP,ublas::matrix<double> initR,ublas::matrix<double> initQ)
{
	cout << "EKF Constructor 2 " << endl;

	n = initX.size1();

	xhat.resize(n,1,false);
	xhat.clear();

	xhat = initX;

	P.resize(n,n,false);
	P.clear();

	P = initP;

	counter = 0;

	//uPrev.resize(n,1,false);
	//uPrev.clear();

	//Q.resize(n,n,false);
	//Q.clear();

	//R.resize(n,n,false);
	//R.clear();

	//	initTuningStates();

	//	timeLinearity = true;//assume linear time update

	//	measLinearity = true;//assume linear measurement update
}

/////////////////////////////////
/*
initStates
initX
- initX is n x 1 column matrix
- each index is filled with initial value
*/
/////////////////////////////////
//void EKF::initStates(ublas::matrix <double> initX)
//{
//	numStates = initX.size1();
//
//	xhat.resize(n,1,false);
//	xhat.clear();
//
//	xhat = initX;
//
//}

void EKF::setP(ublas::matrix<double> Pp)
{
	P.resize(Pp.size1(),Pp.size2());
	P.clear();
	P = Pp;
}

void EKF::setXhat(ublas::matrix<double> xhatp)
{
	xhat.resize(xhatp.size1(),xhatp.size2());
	xhat.clear();
	xhat = xhatp;

}
ublas::matrix<double> EKF::getP()
{
	return P;
}

ublas::matrix<double> EKF::getXhat()
{
	return xhat;
}

/////////////////////////////////
//Deconstructor
EKF::~EKF()
{
	cout << "EKF Deconstructor" << endl;








}

//////////////////////////////////
/*
reset
- reset filter
*/
//////////////////////////////////
void EKF::reset()
{


}

///////////////////////////////////////
/*
getH
- ublas::matrix<double> meas
- only really need the length of meas
- m x 1
- measStateSel
- n x 1
- contains true or false in index whether xhat(index) has a measurement or not
- Ex. [true true false true] means the first 2 states in the state vector and
the last state in the 4 state vector have measurements
- length of meas must equal the number of trues in measStateSel
- meant for filling in H matrices with only 1's and 0's
- returns m x n H matrix

*/
///////////////////////////////////////
//ublas::matrix<double> EKF::getH(int numMeas,ublas::matrix <bool> measStateSel)
//{
//	ublas::matrix<double> returnVar;
//	returnVar.resize(m,n,false);
//	returnVar.clear();
//
//	int counter = 0;
//	for(int i = 0; i<n; i++)
//	{
//		if(measStateSel(i,0) == true)
//		{
//			returnVar(counter,i) = 1;
//			counter = counter + 1;
//		}
//	}
//
//	return returnVar;  
//}

////////////////////////////////////////
/*
rungeKutta
ublas::matrix<double> u
- inputs
- numU x 1
- must know input order
ublas::matrix<double> x
- state vector
- n x n
double dt
- time step
returns xhat after rungeKutta
*/
////////////////////////////////////////

//ublas::matrix<double> EKF::rungeKuttaTime(ublas::matrix<double> u,ublas::matrix<double> x,ublas::matrix<double>(*f)(ublas::matrix<double> u,ublas::matrix<double> x),double dt)
//{
//
//	ublas::matrix<double> Ka;
//	ublas::matrix<double> Kb;
//	ublas::matrix<double> Kc;
//	ublas::matrix<double> Kd;
//
//	Ka.resize(n,1,false);
//	Ka = f(u,x);
//
//	ublas::matrix<double> midval = ((.5*(u - uPrev))+u);
//	Kb = f(midval,(x+.5*dt*Ka));
//
//	Kc = f(midval,(x+.5*dt*Kb));
//
//	Kd = f(u,(x+dt*Kc));
//
//	x = x+(.16666666666666666)*dt*(Ka+2*Kb+2*Kc+Kd);//xminus
//
//	uPrev = u;
//	return x;
//}

////////////////////////////////////////
/*
rungeKuttaMeas
ublas::matrix<double> u
- inputs
- numU x 1
- must know input order
ublas::matrix<double> x
- state vector
- n x n
double dt
- time step
returns xhat after rungeKutta
*/
////////////////////////////////////////
//ublas::matrix<double> EKF::rungeKuttaMeas(ublas::matrix<double> x,ublas::matrix<double>(*h)(ublas::matrix<double> x))
//{
//	h(x);
//	
//	return x;
//}

////////////////////////////////////////////
/*
inverse
- invert matrix using matrix2 type
- probably can do better 
*/
////////////////////////////////////////////
ublas::matrix<double> EKF::inverse(ublas::matrix<double> matrix)
{
	int row = matrix.size1();
	int column = matrix.size2();

	matrix2<double> invMat(row,column);
	ublas::matrix<double> invMatrix(row,column);

	for (int i = 0; i<row; i++)
	{
		for (int j = 0; j<column; j++)
		{
			invMat.setvalue(i,j,matrix.at_element(i,j));
		}
	}

	try
	{
		invMat.invert();

		for (int i = 0; i<row; i++)
		{
			for (int j = 0; j<column; j++)
			{
				bool success;
				double value = 0;
				invMat.getvalue(i,j,value,success);
				invMatrix.insert_element(i,j,value);
			}
		}

		//bool invertedFlag = InvertMatrix(fTickF,invFTickF);
		//fTickF.invert();
	}
	catch( ... )
	{
		cout << "Exception raised" << endl;
		//exceptionThrown = true;

	}

	return invMatrix;

}

void EKF::init(ublas::matrix<double> initX,ublas::matrix<double> initP)
{
	n = initX.size1();

	xhat.resize(n,1,false);
	xhat.clear();

	xhat = initX;

	P.resize(n,n,false);
	P.clear();

	P = initP;
}


////////////////////////////////////////////
/*
timeUpdate
u: ublas::matrix
- m x 1 (number of inputs column vector)
- must know order of states
dt: time step
- double
returnVar
- double*
- n x n
*/
////////////////////////////////////////////
//void EKF::timeUpdate(ublas::matrix<double> uP, double dt)
//{
//
//	u = uP;
//
//	if(timeLinearity)
//	{//assume linear time update
//		A = getA(u,xhat,dt);
//
//		xhat = prod(A,xhat);
//
//		ublas::matrix<double> AP = prod(A,P);
//		//ublas::matrix<double> QW = prod(Q,W);
//		P = prod(AP,trans(A))+ Q;//prod(QW,Q);
//
//	}
//	else
//	{
//		//continuous update
//		xhat = rungeKuttaTime(u,xhat,dt);
//
//		A = getA(u,xhat,dt);
//		ublas::matrix<double> AP = prod(A,P);
//		//ublas::matrix<double> QW = prod(Q,W);
//		P = prod(AP,trans(A))+ Q;//prod(QW,Q);
//	}
//}

////////////////////////////////////////////
/*
measUpdate
meas: ublas::matrix
- m x 1 (number of measurements column vector)
- must know order of measurements
nonLinear: nonLinear
- choose between H matrix with ones or H matrix specified in .h file
*/
////////////////////////////////////////////
//void EKF::measUpdate(ublas::matrix<double> meas, ublas::matrix <bool> measStateSel, double dt)
//{
//	//int count = 0;
//	//count trues
//	//for (int i = 0; i < measStateSel.size1(); i++)
//	//{
//	//	if(measStateSel(i,0) == true)
//	//		count = count + 1;
//	//}
//	if(measLinearity)
//	{
//		H = getH(meas.size1(),measStateSel);
//
//		ublas::matrix<double> intermediate;
//		ublas::matrix<double> HP = prod(H,P);
//
//		ublas::matrix<double> R2;
//		R2.resize(meas.size1(),meas.size1(),false);
//		R2.clear();
//
//		//R2 - m x m matrix
//		//assume R matrix from initialization is n x n to account for different measurements
//		int count = 0;
//		for (int i = 0; i < measStateSel.size1(); i++)
//		{
//			if(measStateSel(i,0) == true)
//			{
//				R2(count,count) = R(i,i);
//				count = count + 1;
//			}
//		}
//
//		intermediate = prod(HP,trans(H)) + R2;
//
//		ublas::matrix<double> invIntermediate;
//		invIntermediate = inverse(intermediate);
//
//		ublas::matrix<double> PHTick=prod(P,trans(H));
//
//		K = prod(PHTick,invIntermediate);
//
//		ublas::matrix<double> res =	meas-prod(H,xhat);
//
//		ublas::matrix<double> HX = prod(H,xhat);
//
//		ublas::matrix<double> Kres = prod(K,res);
//
//		xhat = xhat+prod(K,res);
//
//		ublas::identity_matrix<double> eye (numStates);
//		ublas::matrix<double> IMinusKH = eye-prod(K,H);
//		P = prod(IMinusKH,P);
//	}
//	else
//	{
//		H = getH(meas);//Measurement Jacobian
//
//		ublas::matrix<double> intermediate;
//		ublas::matrix<double> HP = prod(H,P);
//
//		ublas::matrix<double> R2;
//		R2.resize(meas.size1(),meas.size1(),false);
//		R2.clear();
//
//		//R2 - m x m matrix
//		//assume R matrix from initialization is n x n to account for different measurements
//		int count = 0;
//		for (int i = 0; i < measStateSel.size1(); i++)
//		{
//			if(measStateSel(i,0) == true)
//			{
//				R2(count,count) = R(i,i);
//				count = count + 1;
//			}
//		}
//
//		intermediate = prod(HP,trans(H)) + R2;
//
//		ublas::matrix<double> invIntermediate;
//		invIntermediate = inverse(intermediate);
//
//		ublas::matrix<double> PHTick=prod(P,trans(H));
//
//		K = prod(PHTick,invIntermediate);
//
//		ublas::matrix<double> HFunc = rungeKuttaMeas(xhat,dt);//nonlinearity h(xhat)
//
//		ublas::matrix<double> res =	meas-HFunc;
//
//		//ublas::matrix<double> HX = prod(H,xhat);
//
//		//ublas::matrix<double> Kres = prod(K,res);
//
//		xhat = xhat+prod(K,res);
//
//		ublas::identity_matrix<double> eye (n);
//		ublas::matrix<double> IMinusKH = eye-prod(K,H);
//		P = prod(IMinusKH,P);
//
//	}
//}



void EKF::timeUpdate(boost::function3<ublas::matrix<double>,ublas::matrix<double>, ublas::matrix<double>, double>f, ublas::matrix<double> u, ublas::matrix<double> A, ublas::matrix<double> W, ublas::matrix<double> Q, double dt)
{
	//cout << xhat << endl;

	//project state ahead
	//////////////////////////////////////////
	xhat = f(xhat,u,dt);
	//////////////////////////////////////////

	//project error covariance ahead
	//////////////////////////////////////////
	ublas::matrix<double> AP = prod(A,P);

	ublas::matrix<double> WQ = prod(W,Q);

	P = prod(AP,trans(A))+ prod(WQ,trans(W));
	//////////////////////////////////////////

	counter++; // for counter used in array in matlab - record starts with 1
	if(record)
	{
		recordPXhat();
		recordTimeData(u,A,W,Q,dt);
	}

	//cout << xhat << endl;

}
void EKF::measUpdate(boost::function1<ublas::matrix<double>,ublas::matrix<double> >h, ublas::matrix<double> z, ublas::matrix<double> H, ublas::matrix<double> V, ublas::matrix<double> R)
{

	m = z.size1();

	//Calculate Gain K
	////////////////////////////////////////////////////////
	ublas::matrix<double> K;// n x m gain / blending factor
	K.resize(n,m,false);

	ublas::matrix<double> HP = prod(H,P);
	ublas::matrix<double> VR = prod(V,R);

	ublas::matrix<double> intermediate;
	intermediate = prod(HP,trans(H)) + prod(VR,trans(V));

	ublas::matrix<double> invIntermediate;
	invIntermediate = inverse(intermediate);

	ublas::matrix<double> PHTick = prod(P,trans(H));

	K = prod(PHTick,invIntermediate);
	/////////////////////////////////////////////////////////

	//update estimate with measurement
	/////////////////////////////////////////////////////////
	ublas::matrix<double> HFunc = h(xhat);//nonlinearity h(xhat)

	ublas::matrix<double> res =	z-HFunc;

	xhat = xhat+prod(K,res);
	/////////////////////////////////////////////////////////

	//update error covariance
	/////////////////////////////////////////////////////////
	ublas::identity_matrix<double> eye (n);
	ublas::matrix<double> IMinusKH = eye-prod(K,H);
	P = prod(IMinusKH,P);
	/////////////////////////////////////////////////////////

	counter++; // for counter used in array in matlab - record starts with 1
	if(record)
	{
		recordPXhat();
		recordMeasurementData(z,H,V,R);
	}

}


void EKF::writeMatrixToFile(const char* filename,const char* matrixName,ublas::matrix<double> matrix)
{
	ofstream floatfile(filename,ios::out | ios::app);

	if(floatfile.is_open())
	{
		//write name of matrix
		floatfile << matrixName << ":" << "\n";

		//write rows / columns
		for (int i = 0; i<matrix.size1(); i++)
		{
			for(int j = 0; j< matrix.size2(); j++)
			{
				floatfile << matrix(i,j) << "	";
			}
			floatfile << "\n";
		}
	}

	floatfile << "\n";

	floatfile.close();

	return;
}

void EKF::writeColumnMatrixToFile(const char* filename,const char* matrixName,ublas::matrix<double> matrix)
{
	ofstream floatfile(filename,ios::out | ios::app);

	//floatfile << matrixName << "\n";

	if(floatfile.is_open())
	{
		//write name of matrix
		//floatfile << matrixName << ":" << "\n";

		//write rows / columns
		for (int i = 0; i<matrix.size1(); i++)
		{
			floatfile << matrix(i,0) << "	";
		}
	}

	floatfile << "\n";

	floatfile.close();

	return;

}

///////////////////////////////////
//setOptions
//recordAll - true value writes data - very slow processing
//eraseAll - immediately erases data
void EKF::setOptions(bool recordAll,bool eraseAll)
{

	record = recordAll;

	if(eraseAll == true)
		eraseData();

	return;

}

void EKF::recordPXhat()
{
	//writeColumnMatrixToFile("xhat.txt","xhat",xhat);
	//writeColumnMatrixToFile("P.txt","P",P);

	writeDoubleToFile("xhat.txt",xhat(0,0));
	writeDoubleToFile("P.txt",P(0,0));
}

void EKF::recordTimeData(ublas::matrix<double> u, ublas::matrix<double> A, ublas::matrix<double> W, ublas::matrix<double> Q, double dt)//variables exclusive to time update
{
	//u,A,W,Q,dt
	//writeColumnMatrixToFile("u.txt","u",u);
	//writeMatrixToFile("A.txt","A",A);
	//writeMatrixToFile("W.txt","W",W);
	//writeMatrixToFile("Q.txt","Q",Q);
	//writeDoubleToFile("dt.txt",dt);
	//writeDoubleToFile("TimeCount.txt",counter);

	writeDoubleToFile("u0.txt",u(0,0));
	writeDoubleToFile("u1.txt",u(1,0));
	writeDoubleToFile("u2.txt",u(2,0));
	writeDoubleToFile("A.txt",A(0,0));
	writeDoubleToFile("W.txt",W(0,0));
	writeDoubleToFile("Q.txt",Q(0,0));
	writeDoubleToFile("dt.txt",dt);
	writeDoubleToFile("TimeCount.txt",counter);
}

void EKF::recordMeasurementData(ublas::matrix<double> z, ublas::matrix<double> H, ublas::matrix<double> V, ublas::matrix<double> R)//variables exclusive to measurement update
{
	//z, H, V, R
	//writeColumnMatrixToFile("z.txt","z",z);
	//writeMatrixToFile("H.txt","H",H);
	//writeMatrixToFile("V.txt","V",V);
	//writeMatrixToFile("R.txt","R",R);
	//writeDoubleToFile("MeasCount.txt",counter);

	writeDoubleToFile("z.txt",z(0,0));
	writeDoubleToFile("H.txt",H(0,0));
	writeDoubleToFile("V.txt",V(0,0));
	writeDoubleToFile("R.txt",R(0,0));
	writeDoubleToFile("MeasCount.txt",counter);


}

void EKF::eraseData()
{
	//erase text files

	//time update
	clearFile("u.txt");

	clearFile("A.txt");
	clearFile("W.txt");
	clearFile("Q.txt");
	clearFile("dt.txt");
	clearFile("TimeCount.txt");

	//measurement update
	clearFile("z.txt");
	clearFile("H.txt");
	clearFile("V.txt");
	clearFile("R.txt");
	clearFile("MeasCount.txt");

	//both updates
	clearFile("xhat.txt");
	clearFile("P.txt");
}

void EKF::clearFile(const char* filename)
{
	ofstream file(filename,ios::out);
	if(file.is_open())
		file << "";

	file.close();
}

void EKF::writeDoubleToFile(const char* filename,double value)
{
	ofstream floatfile(filename,ios::out | ios::app);
	floatfile.precision(12);

	if(floatfile.is_open())
	{
		floatfile << value << "\n";
	}

	floatfile.close();

	return;
}

/////////////////////////////////////
//loadData
//grab data from file
/////////////////////////////////////
ublas::matrix<double> EKF::loadData(const char* filename)//load IMU data into matrix for SIMULATION
{

	ublas::matrix<double> returnVar;


	double rowCounter = 0;
	double columnCounter = 0;
	string line;
	//text IO operations on filename
	//iterate through text file until desired rowNumber is found

	//-----------------------------------
	//count rows and columns of data
	//---------------------------------
	ifstream dataFile;
	if(filename != NULL)
	{
		dataFile.open(filename,ios::in);		
		if(dataFile.is_open() && !dataFile.eof())
		{
			while(!dataFile.eof())
			{
				getline(dataFile,line);
				if(line != "" && line.at(0) != '%')//ignore commented lines 
				{
					rowCounter++;


					if(columnCounter == 0)//only count rows once - not everytime
					{
						//test to see if line begins with % or empty, ignore if it does
						if(line != "" && line.at(0) != '%')
						{

							char * cstr, *p;

							cstr = new char[line.size()+1];
							strcpy(cstr,line.c_str());

							p=strtok(cstr,"	");
							while(p!= NULL)
							{
								columnCounter++;
								p=strtok(NULL,"	");

							}

							delete[] cstr;
						}
					}
					line.clear();


				}
				else
					line.clear();
			}

		}

		//Parse line - Counter (assumes data stays with same format)
		//-------------------------------
		//tester
		//line.append("15 14 13 12 11 10");

		//-------------------------------


		dataFile.close();

		returnVar.resize(rowCounter, columnCounter);//should be huge
		returnVar.clear();


		line.clear();
		//----------------------------
		//put data into returnVar
		//----------------------------

		double currentRow = 0;
		double currentColumn = 0;
		//line;
		//text IO operations on filename
		//iterate through text file
		ifstream loadFile;
		bool commentedLine = true;
		if(filename != NULL)
		{
			loadFile.open(filename,ios::in);

			while(!loadFile.eof())
			{

				if(loadFile.is_open())
				{
					while(commentedLine && !loadFile.eof())
					{
						line.clear();
						getline(loadFile,line);
						if(line != "" && line.at(0) != '%')
						{
							commentedLine = false;

						}
					}
					commentedLine = true;

				}

				//Parse line
				//-------------------------------
				//tester
				//line.append("15 14 13 12 11 10");

				//-------------------------------

				//test to see if line begins with %
				if(line != "" && line.at(0) != '%')
				{

					char * cstr, *p;

					cstr = new char[line.size()+1];
					strcpy(cstr,line.c_str());

					p=strtok(cstr,"	");
					while(p!= NULL)
					{
						if(currentRow < rowCounter && currentColumn < columnCounter)
						{
							returnVar(currentRow,currentColumn) = strtod(p,NULL);

							double tester = returnVar(currentRow,currentColumn);
						}

						p=strtok(NULL,"	");
						currentColumn++;


					}
					currentColumn = 0;

					delete[] cstr;
				}

				currentRow++;
			}
		}

		loadFile.close();
	}

	return returnVar;
}

vector<double> getIMUData(ublas::matrix<double> IMUDataStorage, int IMURowNumber)//get IMU data from IMU or from IMU matrix
{
	vector<double> imuData;
	imuData.clear();//will return all 0's if no data available for frame

	int columns = IMUDataStorage.size2();
	int rows = IMUDataStorage.size1();

	//check to see that time number is within data
	if(rows > IMURowNumber)
	{
		int counter = 0;
		while(counter < columns)
		{
			imuData.push_back(IMUDataStorage(IMURowNumber,counter));

			counter++;
		}


	}
	else//error
	{

	}


	return imuData;





}



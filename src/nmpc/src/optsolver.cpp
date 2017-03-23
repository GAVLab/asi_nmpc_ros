//

#include "optsolver.h"
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <iterator>
//#include <CMOOSProcessConfiguration>

using namespace std;

namespace WarEagle
{

    OptSolver::OptSolver(): x0Dot(7), x0_k2Vec(7), x0_k3Vec(7), x0_k4Vec(7), xNext_cc(7),xNextVec(7), k1(7), k2(7), k3(7), k4(7), costVec(2), uu(4), uNew(4), costVecPat(2), costVecNew(2), uOpt(4), xVecAug(10), Q(10,10), pat(4,9)//, uPat(4,9)
    {
        // m = 544; //kg from prowler_parameters.m
        m = 20; //kg bullshit number for barbie jeep
        hcg = 0.4; //m from prowler_parameters.m
        Ixx = 435; //kg*m^2 from prowler_parameters.m, Probably bullshit
        cfRoll = (666 / 40) * 1000 * 2 * (pow((hcg - 0.39), 2) / M_PI); //This is just a ballpark figure
        k_phiAux = 384.0*180.0 / M_PI; // Nm/rad
        a = 0.7239; //m from prowler_parameters.m
        b = 0.7239; //m from prowler_parameters.m


        a = 0.25; //m for barbie jeep
        b = 0.25; //m for barbie jeep


        iPsi = 3500; //kg*m^2 from prowler_parameters.m, Probably bullshit though
         iPsi = 350.0; //kg*m^2 bullshit number for barbie jeep
        //Magic Formula Parameters
        B = 9.55;
        C = 1.3;
        D = 6920;
        //Iteration number
        itNum = 0; //Make sure this is zero!!

        // This shit is going on in the mission file now
        // Vmax = 4.75;
        // Vmin = 4.5;

        //Cost function setup
        alpha = 1.0; //I'm assuming these have to be doubles since Q is a double matrix
        beta = 1.0;
        gamma = 1000.0;
        eta = 300.0;
        epsilon = 10.0;

        rwCnt = 10;
        colCnt = 10;

        // delThresh = 37.5*(M_PI/180); //Max steer angle threshold of 37.5 deg for Prowler. May want to check this number with Lowell.
        delThresh = 18*(M_PI/180); //Max steer angle threshold for barbie jeep
        delta_uThresh = 12.5*(M_PI/180); //maximum change in steer angle between iterations

        theta = 0.5;
        lambda = 1;
        del_uThresh = 0.01;
        rollThresh = M_PI / 8;
        rAvoid = 1.5;

        n = 4;
        p = 9;

        xObVec.reserve(100);
        yObVec.reserve(100);
        distObVec.reserve(100);

        obDetected = 0;

        dt = 0.01;


        // x0Dot(7);
        // x0_k2Vec.reserve(7);
        // x0_k3Vec.reserve(7);
        // x0_k4Vec.reserve(7);

        // k1.reserve(7);
        // k2.reserve(7);
        // k3.reserve(7);
        // k4.reserve(7);
        // xNext_cc.reserve(7);
        // xNextVec.reserve(7);

        // costVec.reserve(2);
        // costVecNew.reserve(2);
        // costVecPat.reserve(2);

        // uu.reserve(4);
        // uOpt.reserve(4);


        //Initialize vectors here
        // for (int i;i<7;i++)
        // {
        //     x0Dot[i].push_back(0);
        // }





        //Target and obstacle location info
        //vector<double> target(2); //I will need to store this info elsewhere later...or will I?
        //vector<double> obstacle(2);
        //tLook = 2;

    }

    OptSolver::~OptSolver()
    {
        // Move along, nothing to see here
    }
    vector <double> OptSolver::nonlinVehMod_rk4(double t, vector <double>& x0, int k, vector <double>& u, double tLook)
    {

        //Vehicle Parameters


        Vx = u[3];
        //cout << Vx << "\n";

        //Determine Inputs
        tt = t;
        while (tt>tLook){
            tt = tt - tLook;
        }
        deltaF = (u[0] + (u[1] * (tt / tLook)) + (u[2] * (pow((tt / tLook), 2))));
        dF1 = deltaF;
        dF2 = (u[0] + (u[1] * ((tt + dt) / tLook)) + (u[2] * (pow(((tt + dt) / tLook), 2))));

        uIn[0] = dF1;
        uIn[1] = dF2;

        // dF;
        if (k == 1){
            dF = uIn[0];
        }
        else if (k == 2){
            dF = 0.5*(uIn[0] + uIn[1]);
        }
        else if (k == 3){
            dF = 0.5*(uIn[0] + uIn[1]);
        }
        else if (k == 4){
            dF = uIn[1];
        }


        //Extract Vehicle States
        x = x0[0]; //unused
        y = x0[1]; //unused
        phi = x0[2];
        phiDot = x0[3];
        psi = x0[4];
        psiDot = x0[5];
        Vy = x0[6];

        //Front and rear sideslip angles
        alphaF = dF - atan((Vy + (a*psiDot)) / Vx); //front
        alphaR = atan((-Vy + (b*psiDot)) / Vx); //rear
        //Front and rear tire forces
        Fyf = D*sin(C*atan(B*alphaF)); //front
        Fyr = D*sin(C*atan(B*alphaR)); //rear

        //Vehicle dynamic equations
        psiDotDot = ((a * Fyf * cos(dF)) - (b * Fyr)) / iPsi;
        VyDot = ((Fyf * cos(dF) + Fyr) / m) - (Vx * psiDot);
        betaDot = VyDot / Vx; //I think this is somewhat linearized...
        yDot = (sin(psi)*Vx) + (cos(psi)*Vy);
        xDot = (cos(psi)*Vx) - (sin(psi)*Vy);
        phiDotDot = -((cfRoll / Ixx)*phiDot) - ((k_phiAux / Ixx)*phi) + ((m*Vx*hcg / Ixx)*betaDot) + ((m*Vx*hcg / Ixx)*psiDot);

        //cout <<"xDot ="<< Vx << "\n\n";
        x0DotArr[0] = xDot;
        x0DotArr[1] = yDot;
        x0DotArr[2] = phiDot;
        x0DotArr[3] = phiDotDot;
        x0DotArr[4] = psiDot;
        x0DotArr[5] = psiDotDot;
        x0DotArr[6] = VyDot;
        // x0DotArr[7] = { xDot, yDot, phiDot, phiDotDot, psiDot, psiDotDot, VyDot };
        // for (int i = 0; i<7; i++)
        // {
        //     //cout << x0DotArr[i] << "\n"; //this isn't working
        // }

        // vector <double> x0Dot(7);
        for (int i=0;i<7;i++)
        {
            x0Dot[i] = 0;
        }



        for (int i = 0; i<7; i++){
            fill(x0Dot.begin() + i, x0Dot.begin() + i + 1, x0DotArr[i]);
            //cout << x0Dot[i] << "\n"; //this isn't working
        }

        /*x0Dot[0] = xDot;
        x0Dot[1] = yDot;
        x0Dot[2] = phiDot;
        x0Dot[3] = phiDotDot;
        x0Dot[4] = psiDot;
        x0Dot[5] = psiDotDot;
        x0Dot[6] = VyDot;*/
        return x0Dot;
    }

    vector <double> OptSolver::rk4(vector <double>& x0, double t, vector <double>& u, double tLook){
        itNum++;

        //k1 calculation
        k = 1;
        //double *k1;
        // vector <double> k1(7);
        k1 = OptSolver::nonlinVehMod_rk4(t, x0, k, u, tLook);

        // for(int i;i<7;i++)
        // {
        //     cout << "k1:" << k1[i] << endl;
        // }
        

        //k2 calculation
        k = 2;
        x1 = x0[0] + (0.5*dt*k1[0]);
        x2 = x0[1] + (0.5*dt*k1[1]);
        x3 = x0[2] + (0.5*dt*k1[2]);
        x4 = x0[3] + (0.5*dt*k1[3]);
        x5 = x0[4] + (0.5*dt*k1[4]);
        x6 = x0[5] + (0.5*dt*k1[5]);
        x7 = x0[6] + (0.5*dt*k1[6]);
        x0_k2[0] = x1;
        x0_k2[1] = x2;
        x0_k2[2] = x3;
        x0_k2[3] = x4;
        x0_k2[4] = x5;
        x0_k2[5] = x6;
        x0_k2[6] = x7;
        // x0_k2 = { x1, x2, x3, x4, x5, x6, x7 };
        // vector <double> x0_k2Vec(7);
        for (int i = 0; i<7; i++)
        {
            fill(x0_k2Vec.begin() + i, x0_k2Vec.begin() + i + 1, x0_k2[i]);
            //cout << x0[i] << "\n";
        }
        //double *k2;
        // vector <double> k2(7);
        k2 = OptSolver::nonlinVehMod_rk4(t, x0_k2Vec, k, u, tLook);

        //k3 calculation
        k = 3;
        x1 = x0[0] + (0.5*dt*k2[0]);
        x2 = x0[1] + (0.5*dt*k2[1]);
        x3 = x0[2] + (0.5*dt*k2[2]);
        x4 = x0[3] + (0.5*dt*k2[3]);
        x5 = x0[4] + (0.5*dt*k2[4]);
        x6 = x0[5] + (0.5*dt*k2[5]);
        x7 = x0[6] + (0.5*dt*k2[6]);
        x0_k3[0] = x1;
        x0_k3[1] = x2;
        x0_k3[2] = x3;
        x0_k3[3] = x4;
        x0_k3[4] = x5;
        x0_k3[5] = x6;
        x0_k3[6] = x7;
        // x0_k3 = { x1, x2, x3, x4, x5, x6, x7 };
        // vector <double> x0_k3Vec(7);
        for (int j = 0; j<7; j++){
            fill(x0_k3Vec.begin() + j, x0_k3Vec.begin() + j + 1, x0_k3[j]);
            //cout << x0[i] << "\n";
        }
        //double *k2;
        // vector <double> k3(7);
        k3 = OptSolver::nonlinVehMod_rk4(t, x0_k3Vec, k, u, tLook);
        //double *k3;

        //k4 calculation
        k = 4; //value for k4 calculation
        x1 = x0[0] + dt*k3[0];
        x2 = x0[1] + dt*k3[1];
        x3 = x0[2] + dt*k3[2];
        x4 = x0[3] + dt*k3[3];
        x5 = x0[4] + dt*k3[4];
        x6 = x0[5] + dt*k3[5];
        x7 = x0[6] + dt*k3[6];
        x0_k4[0] = x1;
        x0_k4[1] = x2;
        x0_k4[2] = x3;
        x0_k4[3] = x4;
        x0_k4[4] = x5;
        x0_k4[5] = x6;
        x0_k4[6] = x7;
        // x0_k4 = { x1, x2, x3, x4, x5, x6, x7 };
        // vector <double> x0_k4Vec(7);
        for (int i = 0; i<7; i++){
            fill(x0_k4Vec.begin() + i, x0_k4Vec.begin() + i + 1, x0_k4[i]);
            //cout << x0[i] << "\n";
        }
        //double *k2;
        // vector <double> k4(7);
        k4 = OptSolver::nonlinVehMod_rk4(t, x0_k4Vec, k, u, tLook);
        //double *k4;

    
        for (int i = 0; i<7; i++)
        {
            kGroup[i] = k1[i] + (2 * k2[i]) + (2 * k3[i]) + k4[i];
            xNext[i] = x0[i] + ((dt / 6)*kGroup[i]);
        }
        // vector <double> xNextVec(7);
        for (int i = 0; i<7; i++){
            fill(xNextVec.begin() + i, xNextVec.begin() + i + 1, xNext[i]);
            // cout << x0[i] << "\n";
        }
        return xNextVec;


    }


    vector <double> OptSolver::costCalc(vector <double>& x0, vector <double>& x0_begin, vector <double>& u, double tLook, double rollThresh, double rAvoid, vector <double>& target, vector <double>& obstacle, double uLast, double Vmax, double Vmin)
    {

       // itNum++;

        //This is where terrain info is generated. For actual implementation of this algorithm, this information will be gathered from LiDAR data.
        // double square[2] = { 5., 5. };
        // double terrStart[2] = { 1., -35. };
        // int nlpRow = 9;
        // int nlpCol = 14;
        // MatrixXd roughInfo(nlpRow, nlpCol);
        // MatrixXd gradeInfo(nlpRow, nlpCol);
        // for (int ii = 0; ii<nlpRow; ii++){
        //     for (int kk = 0; kk<nlpCol; kk++){
        //         if (ii == 2 && kk == 4){
        //             roughInfo(ii, kk) = 15;
        //         }
        //         else if (ii == 2 && kk == 5){
        //             roughInfo(ii, kk) = 15;
        //         }
        //         else if (ii == 2 && kk == 8){
        //             roughInfo(ii, kk) = 15;
        //         }
        //         else if (ii == 3 && kk == 4){
        //             roughInfo(ii, kk) = 15;
        //         }
        //         else if (ii == 3 && kk == 5){
        //             roughInfo(ii, kk) = 15;
        //         }
        //         else {
        //             roughInfo(ii, kk) = 2;
        //         }
        //     }

        // }

        // for (int iii = 0; iii<nlpRow; iii++){
        //     for (int kkk = 0; kkk<nlpCol; kkk++){
        //         if (iii == 2 && kkk == 6){
        //             gradeInfo(iii, kkk) = 26.5;
        //         }
        //         else if (iii == 2 && kkk == 7){
        //             gradeInfo(iii, kkk) = 26.5;
        //         }
        //         else {
        //             gradeInfo(iii, kkk) = 0;
        //         }
        //     }

        // }






        // MatrixXd Q(10, 10); //put this in the header
        for (int i = 0; i<rwCnt; i++){
            for (int k = 0; k<colCnt; k++){
                if (i == 2 && k == 2){
                    Q(i, k) = alpha;
                }
                else if (i == 3 && k == 3){
                    Q(i, k) = beta;
                }
                else if (i == 7 && k == 7){
                    Q(i, k) = gamma;
                }
                else if (i == 8 && k == 8){
                    Q(i, k) = eta;
                }
                else if (i == 9 && k == 9){
                    Q(i, k) = epsilon;
                }
                else {
                    Q(i, k) = 0;
                }
            }

        }


        t = 0.0;
        Jtot = 0;
        // vector <double> xNext_cc(7);
        // VectorXd xVecAug(10);


        //Compute Cost
        while (t <= tLook)
        {

            if (t == 0)
            {
                for (int k = 0; k<7; k++)
                {
                    x0[k] = x0_begin[k];
                }
            }

            deltaF = u[0] + (u[1] * (t / tLook)) + (u[2] * (pow((t / tLook), 2)));
            xNext_cc = OptSolver::rk4(x0, t, u, tLook);

            //Determine terrain roughness and grade from vehicle location. I don't need to determine roughness and grade in this manner anymore.
            // int xTerrLoc;
            // int yTerrLoc;

            // if (xNext_cc[0] <= 1){
            //     rough = 0;
            //     grade = 0;
            // }
            // else {
            //     xTerrLoc = ceil((xNext_cc[0] - terrStart[0]) / square[0]);
            //         // cout << "Before terrSTart Access" << endl;

            //     yTerrLoc = ceil((-terrStart[1] - xNext_cc[1]) / square[1]);
            //         // cout << "Before roughInfo" << endl;

            //     rough = roughInfo(xTerrLoc, yTerrLoc);
            //         // cout << "Before gradInfo" << endl;

            //     grade = gradeInfo(xTerrLoc, yTerrLoc);
            //     if (rough >= 10){
            //         Vmax = 2;
            //     }
            // }

            //Hard coding this as zeros for now. Will get this info from LiDAR scans later
            // rough = 0;
            // grade = 0;



            xprime = target[0] - xNext_cc[0];
            yprime = target[1] - xNext_cc[1];
            distTarg = pow((pow(xprime, 2) + pow(yprime, 2)), 0.5); //distance to target

            //Obstacle Detection Stuff
            int g = 0;
            xObVec.clear();
            yObVec.clear();
            distObVec.clear();
            // cout << "obstacle size: " << obstacle.size() << endl;
            for(int i=0;i<obstacle.size();i=i+2)
            {
                xObVec.push_back(obstacle[i]);
                // xObVec[g] = obstacle[i];
                yObVec.push_back(obstacle[i+1]);
                // yObVec[g] = obstacle[i+1];
                // distObVec[g] = pow((pow(xObVec[g], 2) + pow(yObVec[g], 2)), 0.5); //distance to obstacle
                distObVec.push_back(pow((pow(xObVec[g], 2) + pow(yObVec[g], 2)), 0.5));
                // cout << "distObVec " << g << " " << distObVec[g] << endl; 
                g = g+1;
                obDetected = 1;
            }

            if(obDetected == 1)
            {
                minDistOb = *min_element(distObVec.begin(),distObVec.end());
                // cout << "min distance: " << minDistOb << endl;
            }

            

            xOb = obstacle[0] - xNext_cc[0];
            yOb = obstacle[1] - xNext_cc[1];
            distOb = pow((pow(xOb, 2) + pow(yOb, 2)), 0.5); //distance to obstacle

            // cout << "xOb: " << xOb << endl;

            for (int j = 0; j<10; j++)
            {
                if (j<7)
                {
                    xVecAug(j) = xNext_cc[j];
                }
                else if (j == 7)
                {
                    xVecAug(j) = distTarg;
                }
                else if (j == 8)
                {
                    xVecAug(j) = rough;
                }
                else if (j == 9)
                {
                    xVecAug(j) = grade;
                }


            }

            J = 0.5*(xVecAug.transpose())*Q*xVecAug;
            Jtot = Jtot + J;

            // double delThresh = 37.5*(M_PI/180); //Max steer angle threshold of 37.5 deg for Prowler. May want to check this number with Lowell.
            // double delta_uThresh = 5*(M_PI/180); //Maximum allowable change in steer angle set to 20 deg.


            //Constraint stuff
            if (abs(xNext_cc[2]) >= rollThresh)
            {
                conVio = 1;
                t = tLook;
            }
            else if (abs(deltaF) > (delThresh))
            {
                conVio = 1;
                t = tLook;
            }
            else if (distOb < rAvoid)
            {
                conVio = 1;
                t = tLook;
            }
            else if (minDistOb < rAvoid && obDetected == 1) // add this flag for obstacle detection
            {
                conVio = 1;
                t = tLook;
            }
            else if (u[3] > Vmax)
            {
                conVio = 1;
                t = tLook;
            }
            else if (u[3] < Vmin)
            {
                conVio = 1;
                t = tLook;
            }
            else if (grade > 40) //super arbitrary threshold value
            { 
                conVio = 1;
                t = tLook;
            }
            else if (rough > 40){ //super arbitrary threshold value
                conVio = 1;
                t = tLook;
            }
            else if (abs(u[0]-uLast)>delta_uThresh)
            {
                conVio = 1;
                t = tLook;
            }
            else {
                conVio = 0;
            }

            t = t + dt;
            for (int i = 0; i<7; i++){
                fill(x0.begin() + i, x0.begin() + i + 1, xNext_cc[i]);
                //cout << x0[i] << "\n";
            }

        }

        // vector <double> costVec(2);

        fill(costVec.begin(), costVec.begin() + 1, Jtot);
        fill(costVec.begin() + 1, costVec.begin() + 2, conVio);

        obDetected = 0; //reset catch for obstacle detection


        return costVec;
    }


    vector <double> OptSolver::patSearchEick(vector <double>& x0, vector <double>& x0_begin, vector <double>& u, double tLook, vector <double>& target, vector <double>& obstacle, double uLast, double Vmax, double Vmin)
    {

        //Constraint Info
        // double rollThresh = M_PI / 8;
        // double rAvoid = 1.5;

        // int n = 4;
        // int p = 9;
    // cout << "patSearch Eick start." << endl;

        // MatrixXd pat(n, p); //this most of all
        pat <<  1, 0, 0, 0, -1, 0, 0, 0, 0,
                0, 1, 0, 0, 0, -1, 0, 0, 0,
                0, 0, 1, 0, 0, 0, -1, 0, 0,
                0, 0, 0, 1, 0, 0, 0, -1, 0;

        //cout << pat(2,5) << "\n";

        //cout << "pat = " << pat << "\n";

        // double theta = 0.5;
        // int lambda = 1;
        del_u = 1;
        // double del_uThresh = 0.01;
        // double costLow; //I really shouldn't need this, but let's try implementing this
        success = 0;
        //int successPat = 0; //not sure if I need this...
        int k;

        // MatrixXd uPat(n,p);
        // vector <double> uu(4);
        // vector <double> uNew(4);
        // vector <double> costVecNew(2);
        // vector <double> costVecPat(2);

        for (int q = 0; q<7; q++){
            x0_begin[q] = x0[q];
        }


        while (abs(del_u) > del_uThresh){


            costVecNew = OptSolver::costCalc(x0, x0_begin, u, tLook, rollThresh, rAvoid, target, obstacle, uLast, Vmax, Vmin);
            costLow = costVecNew[0];
            k = 0;
            while (k < p){
                for (int j = 0; j<n; j++)
                {
                    fill(uu.begin() + j, uu.begin() + j + 1, u[j] + (del_u*pat(j, k)));
                    //cout << pat(j, k) << "\n";
                    //uu(j) = u(j)+(del_u*pat(j,k)); //this is likely not correct
                }
                // cout << "costVecPat[0]_bef:" << costVecPat[0] << endl;
                costVecPat = OptSolver::costCalc(x0, x0_begin, uu, tLook, rollThresh, rAvoid, target, obstacle, uLast, Vmax, Vmin); //don't know if I can keep uu in this form...
                // cout << "costVecPat[0]_aft:" << costVecPat[0] << endl;

                // I am having a problem storing the lowered cost value here
                if ((costVecPat[0] < costLow) && (costVecPat[1] == 0)){

                    //cout << costVecNew[0] << "\n";
                    costLow = costVecPat[0];

                    // cout << "uNew:" << endl;
                    for (int b = 0; b<4; b++){
                        //fill (uNew.begin()+b,uNew.begin()+b+1,uu[b]);

                        uNew[b] = uu[b];
                        // cout << uNew[b] << endl;

                    }

                    //uNew = uu; // check this

                    if (k <= 6){
                        k = 7;
                    }
                    //cout << k << "\n";
                    success = 1;
                    //successPat = 1;
                }


                k = k + 1;

            }


            //Update for successful/unsuccessful iteration
            if (success == 1){
                for (int jj = 0; jj<n; jj++){
                    //fill (u.begin()+jj,u.begin()+jj+1,uNew[jj]);
                    u[jj] = uNew[jj];

                }
                //cout << u[0] << "\n";
                del_u = lambda*del_u;

            }
            else {
                del_u = theta*del_u;
            }

            //Here is where the catchMat stuff would go

            success = 0;

        }

        // vector <double> uOpt(4);
        for (int jj = 0; jj<n; jj++){
            uOpt[jj] = u[jj];
            // cout << uOpt[jj] << endl;
        }

        return uOpt;

    }





}



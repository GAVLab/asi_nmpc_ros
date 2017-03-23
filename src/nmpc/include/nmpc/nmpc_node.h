#ifndef nmpc_node_H
#define nmpc_node_H

#include <ros/console.h>
#include <ros/ros.h>

#include <string>
#include "optsolver.h"
#include "nmpcStructs.h"

#include <sensor_msgs/Imu.h>
#include <ackermann_msgs/AckermannDrive.h>
#include "ublox_msgs/NavVELNED.h"
#include "ublox_msgs/NavPOSLLH.h"
#include "nmpc/uOpt.h"

#include <eigen3/Eigen/Dense> //seems to be fixed...
using Eigen::MatrixXd;



class nmpc_node  
{
public:
	nmpc_node();
	~nmpc_node();

  

    // vector <double> obsDetected;

    vector <double> obDetectVec;
    
    // Stuff from launch file
    double tLook; //control horizon (s)
    // Minimum and maximum velocity values (m/s)
    double Vmax;
    double Vmin;
    // double obFieldViewx;
    // double obFieldViewy;
    double targetLat;
    double targetLong;
    double targetAlt;
    double obstacleLat;
    double obstacleLong;
    double obstacleAlt;









private:

    // Sensor callbacks
    void gpsNavVelNedCallback(const ublox_msgs::NavVELNED::ConstPtr& gps_msg);
    void gpsLLACallback(const ublox_msgs::NavPOSLLH::ConstPtr& gps_msg);
    void imuCallback(const sensor_msgs::Imu::ConstPtr& imu_msg);
    void uOptCallback(const nmpc::uOpt::ConstPtr& uOpt_msg);

    // Other functions
    void resetFlags();
    void checkForInit(nmpcFlags &flags, gpsData &gps_data);
    void publishContVec(std::vector <double> &uOpt);
    void GenerateControlCommand(gpsData &gps_data, imuData &imu_data);
    int convertSteerCommand(double steerAng);

    std::string node_name;
    ros::NodeHandle nh;
    
    // Subscribers
    ros::Subscriber gpsNavVelNed_sub;
    ros::Subscriber gpsLLA_sub;
    ros::Subscriber imu_sub;
    ros::Subscriber lidar_sub;
    ros::Subscriber nmpc_sub; //subscribe to old messages for initial guesses

    // Publishers
    ros::Publisher nmpc_pub; //publish control input vector
    ros::Publisher teensy_pub; //publish controller commands 

    // Stuff from launch file

    // double tLook; //control horizon (s)
    // // Minimum and maximum velocity values (m/s)
    // double Vmax;
    // double Vmin;
    double obFieldViewx;
    double obFieldViewy;
    // double targetLat;
    // double targetLong;
    // double targetAlt;
    // double obstacleLat;
    // double obstacleLong;
    // double obstacleAlt;

    // Stuff that should be here, and isn't garbage from when I didn't know what I was doing
    bool isInit;
    bool isReinit;

    nmpcFlags flags;
    gpsData gps_data;
    imuData imu_data;

    double velInitThresh;
    // Initializer values
    double init_lat;
    double init_long;
    double init_alt;
    double init_course;
    double courseBias;

    std::vector <double> uPrior;
    double enuVeh[3]; //East north up coordinates
    double enuTar[3]; //East north up coordinates
    double enuOb[3]; //East north up coordinates

    VectorXd target;
    VectorXd obstacle;
    VectorXd vehLoc;
    VectorXd vehVel;

    MatrixXd enu2loc;
    MatrixXd enu2loc_veh;

    std::vector <double> x0;
    std::vector <double> x0_begin;

    // These will get fed into the pattern search because they are nice and standard
    std::vector <double> stdTarVec;
    std::vector <double> stdObVec;

    std::vector <double> uGuess; //put this here to give us some wiggle room on what our initial guess is for uOpt
    std::vector <double> uOpt;

    // Time step for lower level control
    double dt_velControl;
    double dt_steerControl;

    double current_time_velCont;
    double current_time_steerCont;
    double new_uOpt_time;
    // Use these to clock the time taken by the NMPC
    double time1;
    double time2;

    bool publishMsg; //flag to publish commands to teensy

    double dist2tarThresh; // threshold for when to stop the nmpc

    WarEagle::OptSolver solver;













    // vector<float> A;
    // A.reserve(20000);

    
    int obNum; //number of detected obstacles in field of view





    //This is for reading in obstacle locations
    // vector <unsigned char> obstacleVector;
    // vector<float> A;
    int sizeBinary;
    double obFieldView[2]; //2D field of view in which obstacles will be looked for. Any obstacles detected outside this field (local frame) will be ignored



    double Vx; //longitudinal velocity






};

#endif

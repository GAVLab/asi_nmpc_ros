#include "nmpc_node.h"
#include <limits.h>
#include <float.h>
#include <math.h>
#include "GpsFunctions.h"
#include <vector>

// using namespace std;


nmpc_node::nmpc_node() : enu2loc(3,3), enu2loc_veh(3,3), target(3), obstacle(3), vehLoc(3), vehVel(3)
{

	// Setup ros things
	ros::NodeHandle nh("~");
    node_name = ros::this_node::getName();

    nh.getParam("tLook", tLook);
    nh.getParam("Vmax", Vmax);
    nh.getParam("Vmin", Vmin);
    nh.getParam("obFieldViewx", obFieldViewx);
    nh.getParam("obFieldViewy", obFieldViewy);
    nh.getParam("targetLat", targetLat);
    nh.getParam("targetLong", targetLong);
    nh.getParam("targetAlt", targetAlt);
    nh.getParam("obstacleLat", obstacleLat);
    nh.getParam("obstacleLong", obstacleLong);
    nh.getParam("obstacleAlt", obstacleAlt);

    // This info will get changed as development progresses
    gpsNavVelNed_sub = nh.subscribe<ublox_msgs::NavVELNED>("/gps/navvelned", 0, &nmpc_node::gpsNavVelNedCallback, this); //use this message velocity and course
    gpsLLA_sub = nh.subscribe<ublox_msgs::NavPOSLLH>("/gps/navposllh", 0, &nmpc_node::gpsLLACallback, this); //use this message position
    imu_sub = nh.subscribe<sensor_msgs::Imu>("/xbow440_node/imu/data", 0, &nmpc_node::imuCallback, this); 
    nmpc_sub = nh.subscribe<nmpc::uOpt>("/asi_nmpc/contVec", 0, &nmpc_node::uOptCallback, this); //maybe change this type and topic later
  

    nmpc_pub = nh.advertise<nmpc::uOpt>("contVec", 10); //maybe change this type and topic later
    teensy_pub = nh.advertise<ackermann_msgs::AckermannDrive>("ackermann_cmd", 10); //publish command for teensy controller


	// Flags for initialization. Initialize once, then reinitalize when moving for a decent course
	isInit = 0; 
	isReinit = 0;

	resetFlags();

	velInitThresh = 0.4; //m/s threshold on velocity required to initialize course measurements
	gps_data.localVel_long = 0.;
	current_time_velCont = 0.;
	current_time_steerCont = 0.;

	flags.init = false;
	flags.reinit = false;

	for(int i=0;i<3;i++)
	{
		uPrior.push_back(0.);
		uGuess.push_back(0.);
		uOpt.push_back(0.);
	}
	uPrior.push_back((Vmax+Vmin)/2); //initialize first guess of velocity to be in between constraints
	uGuess.push_back((Vmax+Vmin)/2); //initialize first guess of velocity to be in between constraints
	uOpt.push_back(0.); //this is the output, so no worries what I set it to

	for(int n=0;n<2;n++)
	{
		stdTarVec.push_back(0.);
		stdObVec.push_back(0.);
	}

	for(int i=0;i<3;i++)
	{
		target[i] = 0.;
		obstacle[i] = 0.;
		vehLoc[i] = 0.;
		vehVel[i] = 0.;
		gps_data.nedVel.push_back(0.);
	}

	for(int j=0;j<3;j++)
	{
		for(int k=0;k<3;k++)
		{
			enu2loc(j,k) = 0.;
		}
	}

	for(int jj=0;jj<7;jj++)
	{
    x0.push_back(0.);
    x0_begin.push_back(0.);
	}

	// Time steps for lower level controllers
	dt_velControl = 0.2; // Let's put this at 5 Hz for now
	dt_steerControl = 0.5; // let's put this at 2 Hz for now

  dist2tarThresh = 0.5; // m. Stop the nmpc when we are this distance away from target

	publishMsg = false; //set this up as default










	// A.reserve(2000);

	// obDetectVec.reserve(round((obFieldView[0]*obFieldView[1])*2)); //should be much larger than needed 

	obDetectVec.reserve(50000); //should be much larger than needed 




}

nmpc_node::~nmpc_node() //This is never getting called on shutdown, which is not ideal
{
	//Move along, nothing to see here
}


void nmpc_node::gpsNavVelNedCallback(const ublox_msgs::NavVELNED::ConstPtr& gps_msg) 
{

//  std::cout << "gpsnavvelned Callback Initiated" << std::endl;

	// Store velocity information and convert to m/s
	gps_data.nedVel[0] = gps_msg->velN*1e-2;
	gps_data.nedVel[1] = gps_msg->velE*1e-2;
	gps_data.nedVel[2] = gps_msg->velD*1e-2;

	for(int k=0;k<3;k++)
	{
		vehVel[k] = gps_data.nedVel[k];
	}

	gps_data.velMag = pow(pow(vehVel[0],2)+pow(vehVel[1],2)+pow(vehVel[2],2),0.5);

	flags.newVel = true;

	gps_data.course = gps_msg->heading*1e-5; //read in and convert to degrees

	flags.newGpsCourse = true;

	if(isReinit) //convert to local frame once reinitialized
	{
		vehVel = enu2loc*vehVel;
		gps_data.localVel_long = vehVel[0];
	}

	std::cout << "veh vel: " << gps_data.localVel_long << std::endl;

  // Adding this just for testing
  // imu_data.yawRate = rand() % 10 - 5;
  // flags.newYawRate = true;




	// Check for initialization
	checkForInit(flags, gps_data);
	



}

void nmpc_node::gpsLLACallback(const ublox_msgs::NavPOSLLH::ConstPtr& gps_msg)
{
	// Store position information and convert to degrees
	gps_data.lat = gps_msg->lat*1e-7;
	gps_data.lon = gps_msg->lon*1e-7;
	gps_data.alt = gps_msg->height*1e-3; //meters

	flags.newGpsPos = true;

	// Check for initialization
	checkForInit(flags, gps_data);

}


void nmpc_node::imuCallback(const sensor_msgs::Imu::ConstPtr& imu_msg)
{
	imu_data.yawRate = imu_msg->angular_velocity.z;
	flags.newYawRate = true;
	// Maybe add a counter here for better initialization. Might not want to just take the first measurement we get
}

void nmpc_node::uOptCallback(const nmpc::uOpt::ConstPtr& uOpt_msg)
{

	// new_uOpt_time = uOpt_msg->header.stamp.toSec();
	new_uOpt_time = ros::Time::now().toSec();
	for(int i=0;i<uOpt_msg->uOpt.size();i++)
	{
		uPrior[i] = uOpt_msg->uOpt[i];
	}

	time1 = ros::Time::now().toSec();

	GenerateControlCommand(gps_data, imu_data);

	time2 = ros::Time::now().toSec();

  // std::cout << "Computation Time: " << time2-time1 << std::endl;


	ackermann_msgs::AckermannDrive ackermann_command;

	ackermann_command.speed = 35; //m/s Let's just test steering first
	publishMsg = true;

	// std::cout << "time diff vel: " << (new_uOpt_time-current_time_velCont) << std::endl;
	// std::cout << "time diff steer: " << (new_uOpt_time-current_time_steerCont) << std::endl;

	if((new_uOpt_time-current_time_velCont)+(time2-time1) > dt_velControl) //implement velocity control
	{
    // Velocity Control!!!
		// ackermann_command.speed = uOpt[3]; //m/s
		// ackermann_command.speed = 35; //m/s Let's just test steering first
		// // ackermann_command.speed = 0; //m/s Let's just test steering first

		// std::cout << "Vel command!" << std::endl;

		// current_time_velCont = new_uOpt_time;
		// publishMsg = true;

	}
	if((new_uOpt_time-current_time_steerCont)+(time2-time1) > dt_steerControl) //implement steer control
	{
    // Steer Control!!!
		// Convert to number from -100 to 100
		int convertedAng = convertSteerCommand(uOpt[0]); //convert steer angle (radians) to integer for steer motor
		ackermann_command.steering_angle = convertedAng; //radians
		// ackermann_command.steering_angle = rand() % 100 - 100; //radians. Just something for testing
		// ackermann_command.steering_angle = 0;
		std::cout << "Steer command received: " << ackermann_command.steering_angle << std::endl;

		current_time_steerCont = new_uOpt_time;
		publishMsg = true;

	}

	if(publishMsg == true)
	{
		teensy_pub.publish(ackermann_command);
		// std::cout << "Message Publishing!" << std::endl;
		publishMsg = false;
	}

	

}

void nmpc_node::resetFlags()
{
	flags.newGpsPos = false;
    flags.newGpsCourse = false;
    flags.newVel = false;
    flags.newYawRate = false;
}

int nmpc_node::convertSteerCommand(double steerAng)
{
	return round((100/(18*M_PI/180)*steerAng)); // this is for the barbie jeep


}

void nmpc_node::checkForInit(nmpcFlags &flags, gpsData &gps_data)
{
	if(!flags.init && flags.newGpsPos && flags.newGpsCourse && flags.newVel && flags.newYawRate)
	{
		flags.init = true;
		resetFlags();
	}
	else if(!flags.reinit && flags.init && flags.newGpsPos && flags.newGpsCourse && flags.newVel && flags.newYawRate && gps_data.velMag > velInitThresh)
	{
		flags.reinit = true;
		resetFlags(); //is this even necessary?
	}

	if(isInit == 0 && flags.init == true) //time to initialize
	{
		init_lat = gps_data.lat;
		init_long = gps_data.lon;
		init_alt = gps_data.alt;
		init_course = gps_data.course;

		isInit = 1;
		std::cout << "Initializing!" << std::endl;

		// std::pause(5);



		std::cout << std::setprecision(12)<< init_lat << " " << init_long << " " << init_alt << std::endl;

		// Command some velocity to get reinitialization going
		ackermann_msgs::AckermannDrive ackermann_command;
		ackermann_command.speed = 35;
		teensy_pub.publish(ackermann_command);

		// sleep(1000);


	}
	else if(isReinit == 0 && flags.reinit == true) //time to reinitialize
	{
		init_lat = gps_data.lat;
		init_long = gps_data.lon;
		init_alt = gps_data.alt;
		init_course = gps_data.course;
		courseBias = 90 - init_course; //don't get this confused with the MOOS code

		isReinit = 1;
		std::cout << "Reinitializing!" << std::endl;
		// Go ahead and put obstacles and target into local frame
		GpsFunctions::wgslla2enu(enuTar,targetLat,targetLong,targetAlt,init_lat,init_long,init_alt);
		GpsFunctions::wgslla2enu(enuOb,obstacleLat,obstacleLong,obstacleAlt,init_lat,init_long,init_alt);

		// Set up the rotation matrix for target and obstacle
		enu2loc(0,0) = cos(init_course*M_PI/180);
		enu2loc(0,1) = sin(init_course*M_PI/180);
		enu2loc(1,0) = -sin(init_course*M_PI/180);
		enu2loc(1,1) = cos(init_course*M_PI/180);
		enu2loc(2,2) = 1;

		// Set up the rotation matrix for vehicle
		enu2loc_veh(0,0) = cos(courseBias*M_PI/180);
		enu2loc_veh(0,1) = sin(courseBias*M_PI/180);
		enu2loc_veh(1,0) = -sin(courseBias*M_PI/180);
		enu2loc_veh(1,1) = cos(courseBias*M_PI/180);
		enu2loc_veh(2,2) = 1;

		// Set up target and obstacle vectors
		for(int i=0;i<3;i++)
		{
			target[i] = enuTar[i];
			obstacle[i] = enuOb[i];
		}

		// Finally, convert to local frame
		target = enu2loc*target;
		obstacle = enu2loc*obstacle;

		


		// Now let's publish a dummy message to get the ball rolling
		publishContVec(uPrior);
	}
}

void nmpc_node::publishContVec(std::vector <double> &uOpt) //maybe add time to this...?
{
	nmpc::uOpt uOpt_msg;
	for(int i=0;i<uOpt.size();i++)
	{
		uOpt_msg.uOpt.push_back(uOpt[i]);
	}

	nmpc_pub.publish(uOpt_msg);
  double poop = 2;
}

void nmpc_node::GenerateControlCommand(gpsData &gps_data, imuData &imu_data) //add lidar data into here eventually
{
	// First things first, lets get the vehicle in the local reference frame
	GpsFunctions::wgslla2enu(enuVeh,gps_data.lat,gps_data.lon,gps_data.alt,init_lat,init_long,init_alt);
	for(int i=0;i<3;i++)
	{
		vehLoc[i] = enuVeh[i];
	}
	vehLoc = enu2loc_veh*vehLoc;

	// Fill up state vector. First with real values
	x0[0] = vehLoc[0];
	x0_begin[0] = vehLoc[0];
	x0[1] = vehLoc[1];
	x0_begin[1] = vehLoc[1];
	x0[4] = gps_data.course - init_course;
	x0_begin[4] = gps_data.course - init_course;

	// Now with bullshit values
	x0[2] = 0.;
	x0_begin[2] = 0.;
	x0[3] = 0.;
	x0_begin[3] = 0.;
	x0[5] = 0.;
	x0_begin[5] = 0.;
	x0[6] = 0.;
	x0_begin[6] = 0.;

	// Fill up obstacle and target vectors. Neglect altitude
	stdTarVec[0] = target[0];
	stdTarVec[1] = target[1];
	stdObVec[0] = obstacle[0];
	stdObVec[1] = obstacle[1];

	// Set up initial guess. I'll do something more with this later. Probably set this equal to uPrior later
	uGuess[0] = 0;
	uGuess[1] = 0;
	uGuess[2] = 0;
	uGuess[3] = (Vmax+Vmin)/2;

  // Check distance from target. If we are close, let's stop
  double dist2tar = pow(pow(x0[0]-stdTarVec[0],2)+pow(x0[1]-stdTarVec[1],2),0.5);

  if(dist2tar > dist2tarThresh)
  {
    uOpt = solver.patSearchEick(x0, x0_begin, uGuess, tLook, stdTarVec, stdObVec, uPrior[0],Vmax,Vmin); //This is for obstacle detection
  }
  else
  {
    std::cout << "Target Reached" << std::endl;
    for(int j = 0;j<4;j++)
    {
      uOpt[j] = 0.;
    }
  }
	publishContVec(uOpt);

}












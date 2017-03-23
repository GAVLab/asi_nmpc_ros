///////////////////////////////////////////////////////////////////////////
//
// Author: Andrew Eick 
//
// Finally getting this code from MOOS put into ROS
//
//////////////////////////    END_GPL    //////////////////////////////////

#include "nmpc_node.h"

int main(int argc ,char * argv[]) {

  sleep(10);

  ros::init(argc,argv,"nmpc_node");
  nmpc_node node;
  ros::spin();
  return 0;

}

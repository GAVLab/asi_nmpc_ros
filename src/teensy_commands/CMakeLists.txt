cmake_minimum_required(VERSION 2.8.3)
project(teensy_commands)


find_package(catkin REQUIRED COMPONENTS
  roscpp
  sensor_msgs
  std_msgs
  ackermann_msgs
  geometry_msgs
)

catkin_package(
  INCLUDE_DIRS 
  LIBRARIES 
  CATKIN_DEPENDS roscpp sensor_msgs std_msgs ackermann_msgs geometry_msgs 
  EXPORTED_TARGETS 
  DEPENDS serial
)


include_directories(
  ${catkin_INCLUDE_DIRS}
)


find_package(serial REQUIRED)
include_directories(${serial_INCLUDE_DIRS})
link_directories(${serial_LIBRARY_DIRS})


## Declare a C++ executable
 add_executable(command_publisher src/command_publisher.cpp)


## Specify libraries to link a library or executable target against
 target_link_libraries(command_publisher ${catkin_LIBRARIES} ${serial_LIBRARIES})

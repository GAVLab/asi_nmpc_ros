#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include <string>
#include <cstdio>

#include "serial/serial.h"
#include <bitset>

// ROS headers
#include <ros/ros.h>

// Control headers
#include <geometry_msgs/PoseStamped.h>
#include <ackermann_msgs/AckermannDrive.h>
#include <sensor_msgs/Joy.h>


using std::string;
using std::exception;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;


// Mode variables for using joystick to switch between Manual and Auto
    int mode = 0;
    int auto_button = 0;
    int manual_button = 0;

// ROS message containers
    // Ackermann msg from controller
        ackermann_msgs::AckermannDrive ackermann_command;

    // Ackermann msg from joystick
        ackermann_msgs::AckermannDrive ackermann_command_joy_pad;

    // Joystick msg
        sensor_msgs::Joy joy_pad;

// Serial Functions
    void my_sleep(unsigned long milliseconds) {
    #ifdef _WIN32
          Sleep(milliseconds); // 100 ms
    #else
          usleep(milliseconds*1000); // 100 ms
    #endif
    }

    void enumerate_ports()
    {
        vector<serial::PortInfo> devices_found = serial::list_ports();

        vector<serial::PortInfo>::iterator iter = devices_found.begin();

        while( iter != devices_found.end() )
        {
            serial::PortInfo device = *iter++;

            printf( "(%s, %s, %s)\n", device.port.c_str(), device.description.c_str(),
         device.hardware_id.c_str() );
        }
    }

    void print_usage()
    {
        cerr << "Usage: test_serial {-e|<serial port address>} ";
        cerr << "<baudrate> [test string]" << endl;
    }

// ROS message callbacks
    void get_ackermann_command_cb(const ackermann_msgs::AckermannDrive::ConstPtr& msg){
        ackermann_command =  *msg;
    }

     void get_ackermann_command_joy_pad_cb(const ackermann_msgs::AckermannDrive::ConstPtr& msg){
        ackermann_command_joy_pad =  *msg;
    }       

    void joy_pad_cb(const sensor_msgs::Joy::ConstPtr& msg){
        joy_pad = *msg;
        auto_button = joy_pad.buttons[6];
        manual_button = joy_pad.buttons[7];
    }

/* Main program: Read an ackermann message from either joystick or other node publishing the command
   and send it to the serial port
*/
int main(int argc, char **argv)
{
    // Initialize ROS node
        ros::init(argc, argv, "command_streaming_node");
        ros::NodeHandle nh;

    // Subscribers and publishers
        // Navigation node subscriber (used in Auto mode)
            ros::Subscriber setpoint_subscriber = nh.subscribe<ackermann_msgs::AckermannDrive>
                    ("/asi_nmpc/ackermann_cmd", 1, get_ackermann_command_cb);

        // Joystick node subscriber (used in Manual mode): This is a python script I wrote for ROS Joy to produce an ackermann message
            ros::Subscriber setpoint_subscriber_joy_pad = nh.subscribe<ackermann_msgs::AckermannDrive>
                    ("ackermann_cmd_joy_pad", 1, get_ackermann_command_joy_pad_cb);

        // Joystick subscriber
            ros::Subscriber joy_pad_subscriber = nh.subscribe<sensor_msgs::Joy>
                    ("joy", 1, joy_pad_cb);

    // Set ROS loop rate
        ros::Rate rate(20.0);   // 20Hz


    // Display default starting mode on screen
        if(mode == 0)
        {
            cout << "SERIAL COMMANDS NODE: Mode = MANUAL" << endl;

        }
        else
        {
            cout << "SERIAL COMMANDS NODE: Mode = AUTO" << endl;

        }

    // Set serial port
      string port("/dev/ttyUSB0");

    // Set the baudrate
      unsigned long baud = 9600;

    // Open port with baudrate and timeout (in milliseconds)
        serial::Serial my_serial(port, baud, serial::Timeout::simpleTimeout(1000));

        cout << "SERIAL COMMAND NODE: Is the serial port open?";
        if(my_serial.isOpen())
            cout << " Yes." << endl;
        else
        {
            cout << " No." << endl;
        }

    // Serial publishing loop: 
        while(ros::ok())
        {
            ros::spinOnce();

            // Check joy_pad buttons for operating mode
                // cout << "Buttons" << joy_pad.buttons.size() << endl;
                if(auto_button == 1)
                {
                    // Display only if mode changed
                    if(mode == 0)
                    {
                        cout << "SERIAL COMMANDS NODE: Mode = AUTO" << endl;

                    }
                    mode = 1;   // auto mode, left trigger
                }

                if(manual_button == 1)
                {
                    if(mode == 1)
                    {
                        cout << "SERIAL COMMANDS NODE: Mode = MANUAL" << endl;

                    }
                    mode = 0;   // manual mode, right trigger
                }

            // Check for brake
                // if (ackermann_command.speed < 0)
                // {
                //     y[2] = (uint8_t)(-1 * ackermann_command.speed);
                // }
                // else
                // {
                //     y[2] = 0;
                // }

            // Check mode and send command from ackermann msg

                if(mode == 0)                                                   // MANUAL
                {
                    uint8_t y[3];
                    y[0] = (uint8_t)ackermann_command_joy_pad.speed;
                    y[1] = (uint8_t)ackermann_command_joy_pad.steering_angle;
                    y[2] = 0xff;

                    size_t num_bytes = 3;
                    my_serial.write(&y[0], num_bytes);
                }
                else                                                            // AUTO
                {
                    uint8_t y[3];
                    y[0] = (uint8_t)ackermann_command.speed;
                    y[1] = (uint8_t)ackermann_command.steering_angle;
                    y[2] = 0xff;

                    size_t num_bytes = 3;
                    my_serial.write(&y[0], num_bytes);
                }


                // cout << "Sent = " << std::bitset<8>(y[0]) << " " << std::bitset<8>(y[1]) << " " << std::bitset<8>(y[2]) << endl;

            // Sleep
                rate.sleep();
        }

}

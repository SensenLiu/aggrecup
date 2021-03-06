cmake_minimum_required(VERSION 2.8.3)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++0x")
project(offb_posctl)
add_compile_options(-std=c++11)
## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  geometry_msgs
  std_msgs
  mavros
  roscpp
  rospy
  message_generation
)
include_directories("/home/sensenliu/psoptproject/psopt-master/include")
include_directories ( "/home/sensenliu/coinbrew_ipopt/Ipopt/src/Interfaces" "/home/sensenliu/coinbrew_ipopt/Ipopt/src/Common" "/home/sensenliu/coinbrew_ipopt/Ipopt/src/LinAlg")

find_package(PkgConfig REQUIRED)
pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)

#message("------------------ " adolc)

## Generate messages in the 'msg' folder
 add_message_files(
   FILES
   controlstate.msg
 )

### Generate services in the 'srv' folder
## add_service_files(
##   FILES
##   Service1.srv
##   Service2.srv
## )
#
### Generate actions in the 'action' folder
## add_action_files(
##   FILES
##   Action1.action
##   Action2.action
## )
#
## Generate added messages and services with any dependencies listed here
 generate_messages(
   DEPENDENCIES
   geometry_msgs
         std_msgs
 )

#################################################
### Declare ROS dynamic reconfigure parameters ##
#################################################
#
### To declare and build dynamic reconfigure parameters within this
### package, follow these steps:
### * In the file package.xml:
###   * add a build_depend and a exec_depend tag for "dynamic_reconfigure"
### * In this file (CMakeLists.txt):
###   * add "dynamic_reconfigure" to
###     find_package(catkin REQUIRED COMPONENTS ...)
###   * uncomment the "generate_dynamic_reconfigure_options" section below
###     and list every .cfg file to be processed
#
### Generate dynamic reconfigure parameters in the 'cfg' folder
## generate_dynamic_reconfigure_options(
##   cfg/DynReconf1.cfg
##   cfg/DynReconf2.cfg
## )
#
####################################
### catkin specific configuration ##
####################################
### The catkin_package macro generates cmake config files for your package
### Declare things to be passed to dependent projects
### INCLUDE_DIRS: uncomment this if your package contains header files
### LIBRARIES: libraries you create in this project that dependent projects also need
### CATKIN_DEPENDS: catkin_packages dependent projects also need
### DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
#  INCLUDE_DIRS /home/sensenliu/catkin_ws/src/acado/include
#  LIBRARIES offb_posctl
#  CATKIN_DEPENDS geometry_msgs mavros roscpp rospy
#  DEPENDS system_lib
        CATKIN_DEPENDS message_runtime
)

############
### Build ##
############
#
### Specify additional locations of header files
### Your package locations should be listed before other locations
##include_directories(
##include
##${catkin_INCLUDE_DIRS})
#
### Declare a C++ library
## add_library(${PROJECT_NAME}
##   src/${PROJECT_NAME}/offb_posctl.cpp
## )
#
### Add cmake target dependencies of the library
### as an example, code may need to be generated before libraries
### either from message generation or dynamic reconfigure
## add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
#
### Declare a C++ executable
### With catkin_make all packages are built within a single CMake context
### The recommended prefix ensures that target names across packages don't collide
## add_executable(${PROJECT_NAME}_node src/offb_posctl_node.cpp)
#
### Rename C++ executable without prefix
### The above recommended prefix causes long target names, the following renames the
### target back to the shorter version for ease of user use
### e.g. "rosrun someones_pkg node" instead of "rosrun someones_pkg someones_pkg_node"
## set_target_properties(${PROJECT_NAME}_node PROPERTIES OUTPUT_NAME node PREFIX "")
#
### Add cmake target dependencies of the executable
### same as for the library above
## add_dependencies(${PROJECT_NAME}_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
#
### Specify libraries to link a library or executable target against
#
##############
### Install ##
##############
#
## all install targets should use catkin DESTINATION variables
## See http://ros.org/doc/api/catkin/html/adv_user_guide/variables.html
#
### Mark executable scripts (Python etc.) for installation
### in contrast to setup.py, you can choose the destination
## install(PROGRAMS
##   scripts/my_python_script
##   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
## )
#
### Mark executables and/or libraries for installation
## install(TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_node
##   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
##   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
##   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
## )
#
### Mark cpp header files for installation
## install(DIRECTORY include/${PROJECT_NAME}/
##   DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
##   FILES_MATCHING PATTERN "*.h"
##   PATTERN ".svn" EXCLUDE
## )
#
### Mark other files for installation (e.g. launch and bag files, etc.)
## install(FILES
##   # myfile1
##   # myfile2
##   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
## )
#
##############
### Testing ##
##############
#
### Add gtest based cpp test target and link libraries
## catkin_add_gtest(${PROJECT_NAME}-test test/test_offb_posctl.cpp)
## if(TARGET ${PROJECT_NAME}-test)
##   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
## endif()
#
### Add folders to be run by python nosetests
## catkin_add_nosetests(test)
## add_executable(offb_posctl_node src/offb_posctl.cpp)
## target_link_libraries(offb_posctl_node ${catkin_LIBRARIES})
## add_dependencies(offb_posctl offb_posctl_generate_messages_cpp)
include_directories(
        /home/sensenliu/catkin_ws/src/acado/include
        ${catkin_INCLUDE_DIRS}
        src/)
#
add_library(thelib
        src/Parameter.h
        src/Parameter.cpp
        src/PID.h
        src/PID.cpp
        src/FILTER.h
        src/FILTER.cpp)
#

#
#
#
add_executable(offb_posctl
        src/offb_posctl.cpp)
target_link_libraries(offb_posctl ${catkin_LIBRARIES} thelib)
target_link_libraries(offb_posctl
        ${catkin_LIBRARIES}
        /home/sensenliu/catkin_ws/src/acado/lib/libacado_toolkit_s.so
        )
link_directories("lib")
link_libraries("/home/sensenliu/catkin_ws/src/acado/lib/libacado_toolkit_s.so")

add_executable(acado_lag_control src/acado_lag_control.cpp )
target_link_libraries(acado_lag_control ${catkin_LIBRARIES} thelib)
target_link_libraries(acado_lag_control
        ${catkin_LIBRARIES}
        /home/sensenliu/catkin_ws/src/acado/lib/libacado_toolkit_s.so
        )

add_executable(readtrajectoryfile src/readtrajectoryfile.cpp )
target_link_libraries(readtrajectoryfile ${catkin_LIBRARIES} thelib)

add_executable(minitime_psopt_control  src/minitime_psopt_control.cpp)
target_link_libraries(minitime_psopt_control PUBLIC PSOPT PRIVATE PkgConfig::ipopt adolc ${catkin_LIBRARIES} thelib
         )
#add_dependencies(offb_posctl offb_posctl_gencpp)


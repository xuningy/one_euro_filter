cmake_minimum_required(VERSION 3.1)
project(one_euro_filter)

set(CMAKE_CXX_STANDARD 14)

find_package(catkin REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS
)

catkin_python_setup()

include_directories(include ${catkin_INCLUDE_DIRS})
link_directories(${catkin_LIBRARY_DIRS})
add_definitions(${catkin_DEFINITIONS})

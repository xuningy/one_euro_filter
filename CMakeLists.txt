cmake_minimum_required(VERSION 2.8.3)
project(one_euro_filter)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

find_package(catkin REQUIRED)

catkin_package()

catkin_python_setup()

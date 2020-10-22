#
# @author Tobias Weber (orcid: 0000-0002-7230-1932)
# @date 15-aug-2020
# @license: see 'LICENSE' file
#

project(geo)
cmake_minimum_required(VERSION 3.0)


find_package(Qt5Core REQUIRED)
find_package(Qt5Widgets REQUIRED)


set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)


set(CMAKE_CXX_STANDARD 20)
add_definitions(-std=c++20)
add_definitions(-Wall -Wextra -Weffc++)


include_directories("${PROJECT_SOURCE_DIR}")


add_executable(geo
	src/geo_gui.cpp src/geo_gui.h
	src/geo_algos.h src/math_algos.h
)


target_link_libraries(geo
	Qt5::Core Qt5::Widgets
	-lqhull_r -lqhullcpp
)

#
# @author Tobias Weber (orcid: 0000-0002-7230-1932)
# @date 15-aug-2020
# @license: see 'LICENSE' file
#

project(geo)
cmake_minimum_required(VERSION 3.0)


find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5Svg REQUIRED)
find_package(Boost REQUIRED)


set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)


set(CMAKE_CXX_STANDARD 20)
add_compile_options(-std=c++20)
add_compile_options(-Wall -Wextra -Weffc++)


add_definitions(${Boost_CXX_FLAGS})

include_directories(
	${PROJECT_SOURCE_DIR}
	${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}/..
	${Qt5Svg_INCLUDE_DIRS}
)


add_executable(voro
	src/voro_gui.cpp src/voro_gui.h
	src/geo_algos.h src/math_algos.h
)

target_link_libraries(voro
	Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Svg
	-lqhull_r -lqhullcpp
)


add_executable(vis
	src/vis_gui.cpp src/vis_gui.h
	src/geo_algos.h src/math_algos.h
)

target_link_libraries(vis
	Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Svg
)


add_executable(lines
	src/lines_gui.cpp src/lines_gui.h
	src/geo_algos.h src/math_algos.h
)

target_link_libraries(lines
	Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Svg
)


add_executable(caster
	src/caster_gui.cpp src/caster_gui.h
)

target_link_libraries(caster
	Qt5::Core Qt5::Gui Qt5::Widgets
	${Boost_LIBRARIES}
)

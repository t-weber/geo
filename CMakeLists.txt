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
find_package(Threads REQUIRED)
#find_package(GMP REQUIRED)


set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)


set(CMAKE_CXX_STANDARD 20)
add_compile_options(-std=c++20)
add_compile_options(-Wall -Wextra -Weffc++)


add_definitions(${Boost_CXX_FLAGS})

include_directories(
	${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/src
	${PROJECT_SOURCE_DIR}/ext
	${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}/..
	${Qt5Svg_INCLUDE_DIRS}
)


add_executable(voro
	src/gui/voro_gui.cpp src/gui/voro_gui.h
	src/geo_algos.h src/math_algos.h
)
target_link_libraries(voro
	Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Svg
	${Boost_LIBRARIES}
	-lqhull_r -lqhullcpp
	#-lqhullstatic_r -lqhullcpp
)


add_executable(hull
	src/gui/hull_gui.cpp src/gui/hull_gui.h
	src/geo_algos.h src/math_algos.h
)
target_link_libraries(hull
	Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Svg
	${Boost_LIBRARIES}
	-lqhull_r -lqhullcpp
	#-lqhullstatic_r -lqhullcpp
)


add_executable(vis
	src/gui/vis_gui.cpp src/gui/vis_gui.h
	src/geo_algos.h src/math_algos.h
)
target_link_libraries(vis
	Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Svg
	${Boost_LIBRARIES}
)


add_executable(lines
	src/gui/lines_gui.cpp src/gui/lines_gui.h
	src/geo_algos.h src/math_algos.h
)
target_link_libraries(lines
	Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Svg
	${Boost_LIBRARIES}
	Threads::Threads
	-lgmp
)


add_executable(caster
	src/gui/caster_gui.cpp src/gui/caster_gui.h
)
target_link_libraries(caster
	Qt5::Core Qt5::Gui Qt5::Widgets
	${Boost_LIBRARIES}
)


add_executable(closestpair src/cli/closestpair.cpp)
target_link_libraries(closestpair ${Boost_LIBRARIES})


add_executable(closestpair-nd src/cli/closestpair-nd.cpp)
target_link_libraries(closestpair-nd ${Boost_LIBRARIES})


add_executable(maxsubvec src/cli/maxsubvec.cpp)
target_link_libraries(maxsubvec ${Boost_LIBRARIES})


add_executable(sort_by_ch src/cli/sort_by_ch.cpp)
target_link_libraries(sort_by_ch ${Boost_LIBRARIES})

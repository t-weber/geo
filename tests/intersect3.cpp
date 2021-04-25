/**
 * intersection tests
 * @author Tobias Weber
 * @date 25-apr-2021
 * @license see 'LICENSE' file
 *
 * References:
 *  * http://www.boost.org/doc/libs/1_76_0/libs/geometry/doc/html/index.html
 *  * https://www.boost.org/doc/libs/1_76_0/libs/geometry/doc/html/geometry/reference/algorithms/buffer/buffer_7_with_strategies.html
 *  * https://github.com/boostorg/geometry/tree/develop/example
 *
 * g++ -std=c++20 -o intersect3 intersect3.cpp
 */

#include <iostream>
#include <fstream>

#include <boost/geometry.hpp>

namespace geo = boost::geometry;
namespace strat = geo::strategy::buffer;

#include "../src/geo_algos.h"

using t_real = double;
using t_vec = m::vec<t_real, std::vector>;

template<class T = t_real> using t_vertex = geo::model::point<T, 2, geo::cs::cartesian>;
template<class T = t_real> using t_lines = geo::model::linestring<t_vertex<T>, std::vector>;
template<class T = t_real> using t_poly = geo::model::polygon<t_vertex<T>, true /*cw*/, false /*closed*/>;
template<class T = t_real> using t_polys = geo::model::multi_polygon<t_poly<T>, std::vector>;
template<class T = t_real> using t_box = geo::model::box<t_vertex<T>>;
template<class T = t_real> using t_ring = geo::model::ring<t_vertex<T>, true /*cw*/, false /*closed*/, std::vector>;

template<class T = t_real> using t_svg = geo::svg_mapper<t_vertex<T>>;


int main()
{
	// boxes
	t_box<t_real> box1;
	box1.min_corner() = t_vertex<t_real>{-1., 1.};
	box1.max_corner() = t_vertex<t_real>{3., 4.};

	t_box<t_real> box2;
	box2.min_corner() = t_vertex<t_real>{0., 2.};
	box2.max_corner() = t_vertex<t_real>{4., 5.};

	// ------------------------------------------------------------------------
	// intersections
	//std::vector<t_vertex<t_real>> inters_box_box;
	//geo::intersection(box1, box2, inters_box_box);
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// custom intersection calculation
	using t_line = std::tuple<t_vec, t_vec, int>;
	std::vector<t_line> lines
	{{
		std::make_tuple(m::create<t_vec>({-1., 1.}), m::create<t_vec>({3., 1.}), 0),
		std::make_tuple(m::create<t_vec>({3., 1.}), m::create<t_vec>({3., 4.}), 0),
		std::make_tuple(m::create<t_vec>({3., 4.}), m::create<t_vec>({-1., 4.}), 0),
		std::make_tuple(m::create<t_vec>({-1., 4.}), m::create<t_vec>({-1., 1.}), 0),

		std::make_tuple(m::create<t_vec>({0., 2.}), m::create<t_vec>({4., 2.}), 1),
		std::make_tuple(m::create<t_vec>({4., 2.}), m::create<t_vec>({4., 5.}), 1),
		std::make_tuple(m::create<t_vec>({4., 5.}), m::create<t_vec>({0., 5.}), 1),
		std::make_tuple(m::create<t_vec>({0., 5.}), m::create<t_vec>({0., 2.}), 1),
	}};

	auto _custom_inters_box_box = g::intersect_sweep<t_vec, t_line>(lines, 1e-6);

	std::vector<t_vec> custom_inters_box_box;
	for(const auto& tup : _custom_inters_box_box)
		custom_inters_box_box.push_back(std::get<2>(tup));
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// print intersection points
	std::cout << "custom box-box intersection points:" << std::endl;
	for(const t_vec& pt : custom_inters_box_box)
	{
		using namespace m_ops;
		std::cout << "\t" << pt << std::endl;
	}
	std::cout << std::endl;

	/*std::cout << "boost::geo box-box intersection points:" << std::endl;
	for(const auto& vert : inters_box_box)
		std::cout << "\t" << geo::get<0>(vert) << "; " << geo::get<1>(vert) << std::endl;
	std::cout << std::endl;*/

	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// svg
	std::ofstream ofstr("tst.svg");
	t_svg<t_real> svg1{ofstr, 100, 100, "width=\"500px\" height=\"500px\""};

	svg1.add(box1);
	svg1.map(box1, "stroke:#000000; stroke-width:1px; fill:none;", 1.);

	svg1.add(box2);
	svg1.map(box2, "stroke:#000000; stroke-width:1px; fill:none;", 1.);


	/*for(const auto& vert : inters_box_box)
	{
		svg1.add(vert);
		svg1.map(vert, "stroke:#0000ff; stroke-width:1px; fill:#0000ff;", 1.);
		//svg1.text(vert, "inters.", "font-family:\'DejaVu Sans\'; font-size:4pt", 2., 2., 8.);
	}*/

	for(const t_vec& vec : custom_inters_box_box)
	{
		t_vertex<t_real> vert{vec[0], vec[1]};
		svg1.add(vert);
		svg1.map(vert, "stroke:#777700; stroke-width:1px; fill:#777700;", 1.);
	}
	// ------------------------------------------------------------------------

	return 0;
}

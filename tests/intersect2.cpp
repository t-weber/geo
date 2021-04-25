/**
 * intersection tests
 * @author Tobias Weber
 * @date 24-apr-2021
 * @license see 'LICENSE' file
 *
 * References:
 *  * http://www.boost.org/doc/libs/1_76_0/libs/geometry/doc/html/index.html
 *  * https://www.boost.org/doc/libs/1_76_0/libs/geometry/doc/html/geometry/reference/algorithms/buffer/buffer_7_with_strategies.html
 *  * https://github.com/boostorg/geometry/tree/develop/example
 *
 * g++ -std=c++20 -o intersect2 intersect2.cpp
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
	// circle
	t_polys<t_real> circle;
	geo::buffer(t_vertex<t_real>{1., 2.}, circle,
		strat::distance_symmetric<t_real>{2.5},
		strat::side_straight{}, strat::join_round{},
		strat::end_round{}, strat::point_circle{});

	t_polys<t_real> circle2;
	geo::buffer(t_vertex<t_real>{3.5, -0.5}, circle2,
		strat::distance_symmetric<t_real>{3.},
		strat::side_straight{}, strat::join_round{},
		strat::end_round{}, strat::point_circle{});


	// line
	t_lines<t_real> line;
	line.emplace_back(t_vertex<t_real>{-3., -2.});
	line.emplace_back(t_vertex<t_real>{5., 8.});


	// box
	t_box<t_real> box;
	box.min_corner() = t_vertex<t_real>{-1., 1.};
	box.max_corner() = t_vertex<t_real>{3., 4.};


	// ------------------------------------------------------------------------
	// intersections
	std::vector<t_vertex<t_real>> inters_circle_circle;
	geo::intersection(circle, circle2, inters_circle_circle);

	std::vector<t_vertex<t_real>> inters_box_circle;
	geo::intersection(box, circle, inters_box_circle);

	std::vector<t_vertex<t_real>> inters_line_circle;
	geo::intersection(line, circle, inters_line_circle);

	//std::vector<t_vertex<t_real>> inters_box_line;
	//geo::intersection(box, line, inters_box_line);
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// custom intersection calculation
	auto custom_inters_circle_circle = m::intersect_circle_circle<t_vec>(
		m::create<t_vec>({1., 2.}), 2.5,
		m::create<t_vec>({3.5, -0.5}), 3.);

	auto custom_inters_line_circle = m::intersect_line_sphere<t_vec>(
		m::create<t_vec>({-3., -2.}),
		m::create<t_vec>({5., 8.}) - m::create<t_vec>({-3., -2.}),
		m::create<t_vec>({1., 2.}), 2.5,
		false, true);

	std::vector<t_vec> box_pts
	{{
		m::create<t_vec>({-1., 1.}),
		m::create<t_vec>({3., 1.}),
		m::create<t_vec>({3., 4.}),
		m::create<t_vec>({-1., 4.}),
	}};

	auto custom_inters_line_box = g::intersect_line_polylines<t_vec>(
		m::create<t_vec>({-3., -2.}), m::create<t_vec>({5., 8.}), box_pts, true);

	auto custom_inters_box_circle = g::intersect_circle_polylines<t_vec>(
			m::create<t_vec>({1., 2.}), 2.5, box_pts, true);
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// print intersection points
	std::cout << "custom circle-circle intersection points:" << std::endl;
	for(const t_vec& pt : custom_inters_circle_circle)
	{
		using namespace m_ops;
		std::cout << "\t" << pt << std::endl;
	}
	std::cout << std::endl;

	std::cout << "boost::geo circle-circle intersection points:" << std::endl;
	for(const auto& vert : inters_circle_circle)
		std::cout << "\t" << geo::get<0>(vert) << "; " << geo::get<1>(vert) << std::endl;
	std::cout << std::endl;

	std::cout << "custom line-circle intersection points:" << std::endl;
	for(const t_vec& pt : custom_inters_line_circle)
	{
		using namespace m_ops;
		std::cout << "\t" << pt << std::endl;
	}
	std::cout << std::endl;

	std::cout << "boost::geo line-circle intersection points:" << std::endl;
	for(const auto& vert : inters_line_circle)
		std::cout << "\t" << geo::get<0>(vert) << "; " << geo::get<1>(vert) << std::endl;
	std::cout << std::endl;

	std::cout << "custom poly-circle intersection points:" << std::endl;
	for(const t_vec& pt : custom_inters_box_circle)
	{
		using namespace m_ops;
		std::cout << "\t" << pt << std::endl;
	}
	std::cout << std::endl;

	std::cout << "boost::geo poly-circle intersection points:" << std::endl;
	for(const auto& vert : inters_box_circle)
		std::cout << "\t" << geo::get<0>(vert) << "; " << geo::get<1>(vert) << std::endl;
	std::cout << std::endl;
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// svg
	std::ofstream ofstr("tst.svg");
	t_svg<t_real> svg1{ofstr, 100, 100, "width=\"500px\" height=\"500px\""};

	svg1.add(line);
	svg1.map(line, "stroke:#000000; stroke-width:1px; fill:none;", 1.);

	svg1.add(box);
	svg1.map(box, "stroke:#000000; stroke-width:1px; fill:none;", 1.);

	svg1.add(circle);
	svg1.map(circle, "stroke:#007700; stroke-width:1px; fill:none;", 1.);

	svg1.add(circle2);
	svg1.map(circle2, "stroke:#007700; stroke-width:1px; fill:none;", 1.);

	for(const auto& vert : inters_box_circle)
	{
		svg1.add(vert);
		svg1.map(vert, "stroke:#0000ff; stroke-width:1px; fill:#0000ff;", 1.);
		//svg1.text(vert, "inters.", "font-family:\'DejaVu Sans\'; font-size:4pt", 2., 2., 8.);
	}

	for(const auto& vert : inters_line_circle)
	{
		svg1.add(vert);
		svg1.map(vert, "stroke:#007700; stroke-width:1px; fill:#007700;", 1.);
		//svg1.text(vert, "inters.", "font-family:\'DejaVu Sans\'; font-size:4pt", 2., 2., 8.);
	}

	for(const auto& vert : inters_circle_circle)
	{
		svg1.add(vert);
		svg1.map(vert, "stroke:#ff0000; stroke-width:1px; fill:#ff0000;", 1.);
		//svg1.text(vert, "inters.", "font-family:\'DejaVu Sans\'; font-size:4pt", 2., 2., 8.);
	}

	for(const t_vec& vec : custom_inters_line_box)
	{
		t_vertex<t_real> vert{vec[0], vec[1]};
		svg1.add(vert);
		svg1.map(vert, "stroke:#777700; stroke-width:1px; fill:#777700;", 1.);
	}
	// ------------------------------------------------------------------------

	return 0;
}

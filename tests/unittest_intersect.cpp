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
 *  * https://www.boost.org/doc/libs/1_76_0/libs/test/doc/html/index.html
 *
 * g++ -std=c++20 -o unittest_intersect unittest_intersect.cpp
 */

#define BOOST_TEST_MODULE test_intersections

#include <iostream>
#include <fstream>
#include <tuple>

#include <boost/geometry.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/type_index.hpp>
namespace geo = boost::geometry;
namespace strat = geo::strategy::buffer;
namespace test = boost::unit_test;
namespace ty = boost::typeindex;

#include "../src/geo_algos.h"

template<class T> using t_vec = m::vec<T, std::vector>;
template<class T> using t_vertex = geo::model::point<T, 2, geo::cs::cartesian>;
template<class T> using t_poly = geo::model::polygon<t_vertex<T>, true /*cw*/, false /*closed*/>;
template<class T> using t_polys = geo::model::multi_polygon<t_poly<T>, std::vector>;


BOOST_AUTO_TEST_CASE_TEMPLATE(test_intersections, t_real, decltype(std::tuple</*float,*/ double, long double>{}))
{
	std::cout << "Testing with " << ty::type_id_with_cvr<t_real>().pretty_name() << " type." << std::endl;

	constexpr const std::size_t NUM_TESTS = 1000;
	const t_real eps = std::sqrt(std::numeric_limits<t_real>::epsilon());
	const t_real cmp_eps = 1e-2;

	t_real rad_min = 0.5;
	t_real rad_max = 5.;
	t_real x_min = -5.;
	t_real x_max = 5.;
	t_real y_min = -5.;
	t_real y_max = 5.;

	for(std::size_t i=0; i<NUM_TESTS; ++i)
	{
		t_real rad1 = g::get_rand<t_real>(rad_min, rad_max);
		t_real rad2 = g::get_rand<t_real>(rad_min, rad_max);
		t_real x1 = g::get_rand<t_real>(x_min, x_max);
		t_real x2 = g::get_rand<t_real>(x_min, x_max);
		t_real y1 = g::get_rand<t_real>(y_min, y_max);
		t_real y2 = g::get_rand<t_real>(y_min, y_max);


		// circles
		t_polys<t_real> circle1;
		geo::buffer(t_vertex<t_real>{x1, y1},
			circle1,
			strat::distance_symmetric<t_real>{rad1},
			strat::side_straight{},
			strat::join_round{128},
			strat::end_round{128},
			strat::point_circle{512});

		t_polys<t_real> circle2;
		geo::buffer(t_vertex<t_real>{x2, y2},
			circle2,
			strat::distance_symmetric<t_real>{rad2},
			strat::side_straight{},
			strat::join_round{128},
			strat::end_round{128},
			strat::point_circle{512});


		// intersections
		std::vector<t_vertex<t_real>> inters_circle_circle;
		geo::intersection(circle1, circle2, inters_circle_circle);


		// custom intersection calculation
		auto custom_inters = m::intersect_circle_circle<t_vec<t_real>>(
			m::create<t_vec<t_real>>({x1, y1}), rad1,
			m::create<t_vec<t_real>>({x2, y2}), rad2,
			eps);


		auto print_circles = [&x1, &x2, &y1, &y2, &rad1, &rad2, &custom_inters, &inters_circle_circle]() -> void
		{
			std::cout << "--------------------------------------------------------------------------------" << std::endl;
			std::cout << "circle 1: mid = (" << x1 << ", " << y1 << "), rad = " << rad1 << std::endl;
			std::cout << "circle 2: mid = (" << x2 << ", " << y2 << "), rad = " << rad2 << std::endl;
			std::cout << std::endl;

			std::cout << "custom circle-circle intersection points:" << std::endl;
			for(const t_vec<t_real>& pt : custom_inters)
			{
				using namespace m_ops;
				std::cout << "\t" << pt << std::endl;
			}
			std::cout << std::endl;

			std::cout << "boost::geo circle-circle intersection points:" << std::endl;
			for(const auto& vert : inters_circle_circle)
			{
				std::cout << "\t" << geo::get<0>(vert) << "; " << geo::get<1>(vert) << std::endl;
			}
			std::cout << "--------------------------------------------------------------------------------" << std::endl;
			std::cout << std::endl;
		};


		bool sizes_equal = (inters_circle_circle.size() == custom_inters.size());
		BOOST_TEST((sizes_equal));
		if(!sizes_equal)
		{
			print_circles();
			continue;
		}

		if(inters_circle_circle.size() == custom_inters.size() && custom_inters.size() == 1)
		{
			bool equals = m::equals<t_real>(custom_inters[0][0], geo::get<0>(inters_circle_circle[0]), cmp_eps) &&
				m::equals<t_real>(custom_inters[0][1], geo::get<1>(inters_circle_circle[0]), cmp_eps);

			BOOST_TEST((equals));
			if(!equals)
				print_circles();
		}
		else if(inters_circle_circle.size() == custom_inters.size() && custom_inters.size() == 2)
		{
			auto pos1 = std::vector<std::tuple<t_real, t_real>>
			{{
				std::make_tuple(custom_inters[0][0], custom_inters[0][1]),
				std::make_tuple(custom_inters[1][0], custom_inters[1][1]),
			}};

			auto pos2 = std::vector<std::tuple<t_real, t_real>>
			{{
				std::make_tuple(geo::get<0>(inters_circle_circle[0]), geo::get<1>(inters_circle_circle[0])),
				std::make_tuple(geo::get<0>(inters_circle_circle[1]), geo::get<1>(inters_circle_circle[1])),
			}};


			/*std::sort(pos1.begin(), pos1.end(), [](const auto& tup1, const auto& tup2) -> bool
			{
				return std::get<0>(tup1) < std::get<0>(tup2);
			});*/

			std::sort(pos2.begin(), pos2.end(), [](const auto& tup1, const auto& tup2) -> bool
			{
				return std::get<0>(tup1) < std::get<0>(tup2);
			});


			bool equals = (m::equals<t_real>(std::get<0>(pos1[0]), std::get<0>(pos2[0]), cmp_eps) &&
				m::equals<t_real>(std::get<1>(pos1[0]), std::get<1>(pos2[0]), cmp_eps) &&
				m::equals<t_real>(std::get<0>(pos1[1]), std::get<0>(pos2[1]), cmp_eps) &&
				m::equals<t_real>(std::get<1>(pos1[1]), std::get<1>(pos2[1]), cmp_eps));

			BOOST_TEST((equals));
			if(!equals)
				print_circles();
		}
	}
}

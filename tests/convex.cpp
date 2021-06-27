/**
 * split a concave polygon into convex regions
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 26-jun-21
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -Wall -Wextra -o convex convex.cpp
 */

#include <boost/geometry.hpp>
namespace geo = boost::geometry;
namespace strat = geo::strategy::buffer;

#include "../src/geo_algos.h"
using namespace m_ops;


using t_real = double;
using t_vec = m::vec<t_real, std::vector>;


template<class T = t_real> using t_vertex = geo::model::point<T, 2, geo::cs::cartesian>;
template<class T = t_real> using t_lines = geo::model::linestring<t_vertex<T>, std::vector>;
template<class T = t_real> using t_svg = geo::svg_mapper<t_vertex<T>>;


void write_svg(const char* filename,
	const std::vector<t_vec>& orig_poly,
	const std::vector<std::vector<t_vec>>& split_poly)
{
	std::ofstream ofstr(filename);
	t_svg<t_real> svg{ofstr, 100, 100, "width=\"500px\" height=\"500px\""};

	// original polygon
	for(std::size_t vertidx1=0; vertidx1<orig_poly.size(); ++vertidx1)
	{
		std::size_t vertidx2 = (vertidx1+1) % orig_poly.size();

		const t_vec& _vert1 = orig_poly[vertidx1];
		const t_vec& _vert2 = orig_poly[vertidx2];

		t_vertex<t_real> vert1{_vert1[0], _vert1[1]};
		t_vertex<t_real> vert2{_vert2[0], _vert2[1]};

		t_lines<t_real> line;
		line.emplace_back(vert1);
		line.emplace_back(vert2);

		svg.add(line);
		svg.map(line, "stroke:#ff0000; stroke-width:3px; fill:none;", 1.);

		svg.add(vert1);
		svg.map(vert1, "stroke:#ff0000; stroke-width:1px; fill:#ff0000;", 3.);
	}

	// split polygons
	for(const auto& poly : split_poly)
	{
		for(std::size_t vertidx1=0; vertidx1<poly.size(); ++vertidx1)
		{
			std::size_t vertidx2 = (vertidx1+1) % poly.size();

			const t_vec& _vert1 = poly[vertidx1];
			const t_vec& _vert2 = poly[vertidx2];

			t_vertex<t_real> vert1{_vert1[0], _vert1[1]};
			t_vertex<t_real> vert2{_vert2[0], _vert2[1]};

			t_lines<t_real> line;
			line.emplace_back(vert1);
			line.emplace_back(vert2);

			svg.add(line);
			svg.map(line, "stroke:#000000; stroke-width:1px; fill:none;", 1.);

			svg.add(vert1);
			svg.map(vert1, "stroke:#0000ff; stroke-width:1px; fill:#0000ff;", 1.);
		}
	}
}


int main()
{
	{
		std::vector<t_vec> poly
		{{
			m::create<t_vec>({0., 0.}),
			m::create<t_vec>({10., 0.}),
			m::create<t_vec>({5., 5.}),
		}};

		std::cout << "poly1" << std::endl;
		auto split = g::convex_split<t_vec>(poly);
		if(split.size() == 0)
			std::cout << "already convex" << std::endl;
		std::cout << std::endl;
	}


	{
		std::vector<t_vec> poly
		{{
			m::create<t_vec>({0., 0.}),
			m::create<t_vec>({10., 0.}),
			m::create<t_vec>({5., 5.}),
			m::create<t_vec>({5., 10.}),
		}};

		std::cout << "poly2" << std::endl;
		auto split = g::convex_split<t_vec>(poly);
		if(split.size() == 0)
			std::cout << "already convex" << std::endl;

		for(std::size_t idx=0; idx<split.size(); ++idx)
		{
			const auto& polysplit = split[idx];
			std::cout << "split polygon " << idx << std::endl;

			for(const auto& vec : polysplit)
				std::cout << vec << std::endl;
		}
		std::cout << std::endl;
	}


	{
		std::vector<t_vec> poly
		{{
			/*0*/ m::create<t_vec>({0., 0.}),
			/*1*/ m::create<t_vec>({2., 2.}),
			/*2*/ m::create<t_vec>({8., 2.}),
			/*3*/ m::create<t_vec>({10., 0.}),
			/*4*/ m::create<t_vec>({10., 10.}),
			/*5*/ m::create<t_vec>({8., 8.}),
			/*6*/ m::create<t_vec>({2., 8.}),
			/*7*/ m::create<t_vec>({0., 10.}),
		}};

		std::cout << "poly3" << std::endl;
		auto split = g::convex_split<t_vec>(poly, 1e-5);
		if(split.size() == 0)
			std::cout << "already convex" << std::endl;

		for(std::size_t idx=0; idx<split.size(); ++idx)
		{
			const auto& polysplit = split[idx];
			std::cout << "split polygon " << idx << std::endl;

			for(const auto& vec : polysplit)
				std::cout << vec << std::endl;
		}
		std::cout << std::endl;
	}


	{
		std::vector<t_vec> poly
		{{
			/*0*/ m::create<t_vec>({5., 8.}),
			/*1*/ m::create<t_vec>({6.7333333333333, 7.7111111111111}),
			/*2*/ m::create<t_vec>({4., 15.}),
			/*3*/ m::create<t_vec>({-1., 9.}),
		}};

		std::cout << "poly4" << std::endl;
		auto split = g::convex_split<t_vec>(poly, 1e-6, true);
		if(split.size() == 0)
			std::cout << "already convex" << std::endl;

		for(std::size_t idx=0; idx<split.size(); ++idx)
		{
			const auto& polysplit = split[idx];
			std::cout << "split polygon " << idx << std::endl;

			for(const auto& vec : polysplit)
				std::cout << vec << std::endl;
		}
		write_svg("poly4.svg", poly, split);
		std::cout << std::endl;
	}


	//if(0)
	{
		std::vector<t_vec> poly
		{{
			m::create<t_vec>({0., 0.}),
			m::create<t_vec>({3.25, -4.1}),
			m::create<t_vec>({10., 1.}),
			m::create<t_vec>({7., 3.}),
			m::create<t_vec>({9.5, 5.}),
			m::create<t_vec>({15., 6.5}),
			m::create<t_vec>({11., 10.}),
			m::create<t_vec>({7., 7.}),
			m::create<t_vec>({4., 15.}),
			m::create<t_vec>({-1., 9.}),
			m::create<t_vec>({5., 8.}),
			m::create<t_vec>({4., 4.}),
		}};

		std::cout << "poly5" << std::endl;
		auto split = g::convex_split<t_vec>(poly, 1e-6, true);
		if(split.size() == 0)
			std::cout << "already convex" << std::endl;

		for(std::size_t idx=0; idx<split.size(); ++idx)
		{
			const auto& polysplit = split[idx];
			std::cout << "split polygon " << idx << std::endl;

			for(const auto& vec : polysplit)
				std::cout << vec << std::endl;
		}
		write_svg("poly5.svg", poly, split);
		std::cout << std::endl;
	}

	return 0;
}

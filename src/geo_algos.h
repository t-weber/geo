/**
 * geometric calculations
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date Oct/Nov-2020
 * @license: see 'LICENSE' file
 *
 * Reference for the algorithms:
 *   - (Klein 2005) "Algorithmische Geometrie" (2005), ISBN: 978-3540209560 (http://dx.doi.org/10.1007/3-540-27619-X).
 *   - (FUH 2020) "Algorithmische Geometrie" (2020), Kurs 1840, Fernuni Hagen (https://vu.fernuni-hagen.de/lvuweb/lvu/app/Kurs/1840).
 *   - (Berg 2008) "Computational Geometry" (2008), ISBN: 978-3-642-09681-5 (http://dx.doi.org/10.1007/978-3-540-77974-2).
 */

#ifndef __GEO2D_ALGOS_H__
#define __GEO2D_ALGOS_H__

#include <vector>
#include <queue>
#include <list>
#include <tuple>
#include <stack>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <ranges>
#include <limits>
#include <random>
#include <iostream>
#include <fstream>

#include "math_algos.h"
#include "graph_algos.h"
#include "math_conts.h"
#include "geo_conts.h"
#include "graph_conts.h"
#include "helpers.h"

#include <boost/intrusive/bstree.hpp>
#include <boost/intrusive/avltree.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#if !__has_include(<voronoi_visual_utils.hpp>)
	#pragma message("Boost.Polygon's voronoi_visual_utils2.hpp utility header was not found, disabling.")
#else
	#include <boost/polygon/polygon.hpp>
	#include <boost/polygon/voronoi.hpp>
	#include <voronoi_visual_utils.hpp>

	#define __GEO2D_USE_BOOST_POLY__
#endif

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullRidge.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullFacetSet.h>
#include <libqhullcpp/QhullVertexSet.h>



// ----------------------------------------------------------------------------
// make point and line segment classes known for boost.polygon
// @see https://www.boost.org/doc/libs/1_75_0/libs/polygon/doc/gtl_custom_point.htm
// @see https://github.com/boostorg/polygon/blob/develop/example/voronoi_basic_tutorial.cpp
// ----------------------------------------------------------------------------
#ifdef __GEO2D_USE_BOOST_POLY__

template<class t_vec> requires m::is_vec<t_vec>
struct boost::polygon::geometry_concept<t_vec>
{
	using type = boost::polygon::point_concept;
};


template<class t_vec> requires m::is_vec<t_vec>
struct boost::polygon::geometry_concept<std::pair<t_vec, t_vec>>
{
	using type = boost::polygon::segment_concept;
};


template<class t_vec> requires m::is_vec<t_vec>
struct boost::polygon::point_traits<t_vec>
{
	using coordinate_type = typename t_vec::value_type;

	static coordinate_type get(const t_vec& vec, boost::polygon::orientation_2d orientation)
	{
		return vec[orientation.to_int()];
	}
};


template<class t_vec> requires m::is_vec<t_vec>
struct boost::polygon::segment_traits<std::pair<t_vec, t_vec>>
{
	using coordinate_type = typename t_vec::value_type;
	using point_type = t_vec;
	using line_type = std::pair<t_vec, t_vec>; // for convenience, not part of interface

	static const point_type& get(const line_type& line, boost::polygon::direction_1d direction)
	{
		switch(direction.to_int())
		{
			case 1: return std::get<1>(line);
			case 0: default: return std::get<0>(line);
		}
	}
};

#endif
// ----------------------------------------------------------------------------



// geo
namespace g{

// ----------------------------------------------------------------------------
// helper functions
// ----------------------------------------------------------------------------

template<class t_num>
t_num get_rand(t_num min=1, t_num max=-1)
{
	static std::mt19937 rng{std::random_device{}()};

	if(max <= min)
	{
		min = std::numeric_limits<t_num>::lowest() / 10.;
		max = std::numeric_limits<t_num>::max() / 10.;
	}

	if constexpr(std::is_integral_v<t_num>)
		return std::uniform_int_distribution<t_num>(min, max)(rng);
	else
		return std::uniform_real_distribution<t_num>(min, max)(rng);
}


template<class t_vec>
t_vec calc_circumcentre(const std::vector<t_vec>& triag)
requires m::is_vec<t_vec>
{
	using namespace m_ops;

	using t_real = typename t_vec::value_type;

	if(triag.size() < 3)
		return t_vec{};

	const t_vec& v0 = triag[0];
	const t_vec& v1 = triag[1];
	const t_vec& v2 = triag[2];

	// formula, see: https://de.wikipedia.org/wiki/Umkreis
	const t_real x =
		(v0[0]*v0[0]+v0[1]*v0[1]) * (v1[1]-v2[1]) +
		(v1[0]*v1[0]+v1[1]*v1[1]) * (v2[1]-v0[1]) +
		(v2[0]*v2[0]+v2[1]*v2[1]) * (v0[1]-v1[1]);

	const t_real y =
		(v0[0]*v0[0]+v0[1]*v0[1]) * (v2[0]-v1[0]) +
		(v1[0]*v1[0]+v1[1]*v1[1]) * (v0[0]-v2[0]) +
		(v2[0]*v2[0]+v2[1]*v2[1]) * (v1[0]-v0[0]);

	const t_real n =
		t_real{2}*v0[0] * (v1[1]-v2[1]) +
		t_real{2}*v1[0] * (v2[1]-v0[1]) +
		t_real{2}*v2[0] * (v0[1]-v1[1]);

	return m::create<t_vec>({x/n, y/n});
}



template<class t_vec, class t_real=typename t_vec::value_type>
t_real line_angle(const t_vec& pt1, const t_vec& pt2)
requires m::is_vec<t_vec>
{
	t_vec dir = pt2 - pt1;
	return std::atan2(dir[1], dir[0]);
}


template<class t_vec, class t_real=typename t_vec::value_type>
t_real line_angle(const t_vec& line1vert1, const t_vec& line1vert2,
	const t_vec& line2vert1, const t_vec& line2vert2)
requires m::is_vec<t_vec>
{
	return line_angle<t_vec, t_real>(line2vert1, line2vert2)
		- line_angle<t_vec, t_real>(line1vert1, line1vert2);
}


template<class t_vec> requires m::is_vec<t_vec>
std::pair<bool, t_vec> intersect_lines(
	const t_vec& pos1a, const t_vec& pos1b,
	const t_vec& pos2a, const t_vec& pos2b,
	bool only_segments = true,
	typename t_vec::value_type eps=std::numeric_limits<typename t_vec::value_type>::epsilon())
{
	t_vec dir1 = pos1b - pos1a;
	t_vec dir2 = pos2b - pos2a;

	auto[pt1, pt2, valid, dist, param1, param2] =
		m::intersect_line_line(pos1a, dir1, pos2a, dir2, eps);

	if(!valid)
		return std::make_pair(false, m::create<t_vec>({}));

	if(only_segments && (param1<0. || param1>=1. || param2<0. || param2>=1.))
		return std::make_pair(false, m::create<t_vec>({}));

	/*std::cout.precision(20);
	std::cout << "intersection between line segment ";
	print_line<t_vec>(std::cout, std::make_tuple(pos1a, pos1b));

	std::cout << " and ";
	print_line<t_vec>(std::cout, std::make_tuple(pos2a, pos2b));

	std::cout << ": ";
	print_point<t_vec>(std::cout, pt1);
	std::cout << ", ";
	print_point<t_vec>(std::cout, pt2);

	std::cout << ", eps = " << eps << ".";
	std::cout << std::endl;*/

	// check if the intersection points on the two lines are the same
	// to rule out numeric instability
	bool alternatives_equal = m::equals<t_vec>(pt1, pt2, eps);
	return std::make_pair(alternatives_equal, pt1);
}


/**
 * intersection of a line and polygon line segments
 */
template<class t_vec, template<class...> class t_cont=std::vector>
t_cont<t_vec> intersect_line_polylines(
	const t_vec& linePt1, const t_vec& linePt2,
	const t_cont<t_vec>& poly, bool only_segment = false,
	typename t_vec::value_type eps=std::numeric_limits<typename t_vec::value_type>::epsilon())
requires m::is_vec<t_vec>
{
	t_cont<t_vec> inters;

	for(std::size_t idx=0; idx<poly.size(); ++idx)
	{
		std::size_t idx2 = (idx+1) % poly.size();
		const t_vec& pt1 = poly[idx];
		const t_vec& pt2 = poly[idx2];

		if(auto [has_inters, inters_pt] =
			intersect_lines<t_vec>(linePt1, linePt2, pt1, pt2, only_segment, eps);
			has_inters)
		{
			inters.emplace_back(std::move(inters_pt));
		}
	}

	// sort intersections by x
	std::sort(inters.begin(), inters.end(), [](const t_vec& vec1, const t_vec& vec2) -> bool
	{
		return vec1[0] < vec2[0];
	});
	return inters;
}


/**
 * intersection of a circle and polygon line segments
 */
template<class t_vec, template<class...> class t_cont=std::vector>
t_cont<t_vec> intersect_circle_polylines(
	const t_vec& circleOrg, typename t_vec::value_type circleRad,
	const t_cont<t_vec>& poly, bool only_segment = false)
requires m::is_vec<t_vec>
{
	t_cont<t_vec> inters;

	for(std::size_t idx=0; idx<poly.size(); ++idx)
	{
		std::size_t idx2 = (idx+1) % poly.size();
		const t_vec& pt1 = poly[idx];
		const t_vec& pt2 = poly[idx2];

		auto theinters = m::intersect_line_sphere<t_vec>(
			pt1, pt2-pt1, circleOrg, circleRad, false, only_segment);

		for(const t_vec& vec : theinters)
			inters.emplace_back(std::move(vec));
	}

	// sort intersections by x
	std::sort(inters.begin(), inters.end(), [](const t_vec& vec1, const t_vec& vec2) -> bool
	{
		return vec1[0] < vec2[0];
	});

	return inters;
}


template<class t_vec, class t_line=std::pair<t_vec, t_vec>, class t_real=typename t_vec::value_type>
requires m::is_vec<t_vec>
std::pair<t_real, t_real> get_line_slope_offs(const t_line& line)
{
	const t_vec& pt1 = std::get<0>(line);
	const t_vec& pt2 = std::get<1>(line);

	if(pt1.size() < 2 || pt2.size() < 2)
		return std::make_pair(0, 0);

	t_real slope = (pt2[1] - pt1[1]) / (pt2[0] - pt1[0]);
	t_real offs = pt1[1] - pt1[0]*slope;

	return std::make_pair(slope, offs);
}


template<class t_vec, class t_line=std::pair<t_vec, t_vec>, class t_real=typename t_vec::value_type>
requires m::is_vec<t_vec>
t_real get_line_y(const t_line& line, t_real x)
{
	auto [slope, offs] = get_line_slope_offs<t_vec>(line);
	return slope*x + offs;
}


/**
 * are two lines equal?
 */
template<class t_vec, class t_line=std::pair<t_vec, t_vec>, class t_real=typename t_vec::value_type>
requires m::is_vec<t_vec>
bool is_line_equal(const t_line& line1, const t_line& line2,
	t_real eps=std::numeric_limits<t_real>::epsilon(),
	bool test_literal_equality = false)
{
	if(test_literal_equality)
	{
		if(!m::equals<t_vec>(std::get<0>(line1), std::get<0>(line2), eps))
			return false;
		if(!m::equals<t_vec>(std::get<1>(line1), std::get<1>(line2), eps))
			return false;
		return true;
	}
	else
	{
		auto [slope1, offs1] = get_line_slope_offs<t_vec>(line1);
		auto [slope2, offs2] = get_line_slope_offs<t_vec>(line2);

		return m::equals<t_real>(slope1, slope2, eps) && m::equals<t_real>(offs1, offs2, eps);
	}
}


template<class t_line, class t_vec=typename std::tuple_element<0, t_line>::type>
std::tuple<bool, t_vec>
intersect_lines(const t_line& line1, const t_line& line2,
	typename t_vec::value_type eps=std::numeric_limits<typename t_vec::value_type>::epsilon())
requires m::is_vec<t_vec>
{
	return intersect_lines<t_vec>(
		std::get<0>(line1), std::get<1>(line1),
		std::get<0>(line2), std::get<1>(line2),
		true, eps);
}


/**
 * returns > 0 if point is on the left-hand side of line
 */
template<class t_vec, class t_real = typename t_vec::value_type>
t_real side_of_line(const t_vec& vec1a, const t_vec& vec1b, const t_vec& pt)
requires m::is_vec<t_vec>
{
	//return line_angle<t_vec>(vec1a, pt) - line_angle<t_vec>(vec1a, vec1b);
	using namespace m_ops;

	t_vec dir1 = vec1b - vec1a;
	t_vec dir2 = pt - vec1a;

	return dir1[0]*dir2[1] - dir1[1]*dir2[0];
}


template<class t_vec, class t_real = typename t_vec::value_type>
bool all_points_on_same_side(const t_vec& lineA, const t_vec& lineB,
	const std::vector<t_vec>& hullvertices, t_real eps=1e-5)
requires m::is_vec<t_vec>
{
	// find a reference vertex which is sufficiently far from the line
	std::optional<t_real> side;
	for(const t_vec& vert : hullvertices)
	{
		if(!side)
		{
			t_real curside = side_of_line(lineA, lineB, vert);
			if(std::abs(curside) > eps)
				side = curside;
		}

		if(side)
			break;
	}

	if(!side)
		return true;


	// are all other vertices on the same side as the reference vertex (or on the line)?
	for(const t_vec& vert : hullvertices)
	{
		t_real curside = side_of_line(lineA, lineB, vert);
		if(std::signbit(*side) != std::signbit(curside) && std::abs(curside) > eps)
			return false;
	}

	return true;
}


template<class t_vec>
bool pt_inside_hull(const std::vector<t_vec>& hull, const t_vec& pt)
requires m::is_vec<t_vec>
{
	//auto [hull, midpt] = sort_vertices_by_angle<t_vec>(_hull);

	// iterate vertices
	for(std::size_t idx1=0; idx1<hull.size(); ++idx1)
	{
		std::size_t idx2 = idx1+1;
		if(idx2 >= hull.size())
			idx2 = 0;

		const t_vec& vert1 = hull[idx1];
		const t_vec& vert2 = hull[idx2];

		// outside?
		if(side_of_line<t_vec>(vert1, vert2, pt) < 0.)
			return false;
	}

	return true;
}


/**
 * get barycentric coordinates of a point
 * see: https://en.wikipedia.org/wiki/Barycentric_coordinate_system
 */
template<class t_vec>
std::optional<t_vec> get_barycentric(const t_vec& tri1, const t_vec& tri2, const t_vec& tri3, const t_vec& pt)
requires m::is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;
	using t_mat = m::mat<t_real, std::vector>;

	t_mat trafo = m::create<t_mat, t_vec>({tri1-tri3, tri2-tri3});
	auto [inv_trafo, ok] = m::inv<t_mat, t_vec>(trafo);
	if(!ok)
		return std::nullopt;

	t_vec vecBary = inv_trafo * (pt-tri3);

	//using namespace m_ops;
	//std::cout << "trafo: " << trafo << "\ninv_trafo: " << inv_trafo << std::endl;
	//std::cout << "bary: " << vecBary << std::endl;
	return vecBary;
}


template<class t_vec>
bool pt_inside_triag(const t_vec& tri1, const t_vec& tri2, const t_vec& tri3, const t_vec& pt)
requires m::is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;
	auto vecBary = get_barycentric<t_vec>(tri1, tri2, tri3, pt);
	if(!vecBary)
		return false;

	t_real x = (*vecBary)[0];
	t_real y = (*vecBary)[1];
	t_real z = t_real{1} - x - y;

	return x>=t_real{0} && x<t_real{1} &&
		y>=t_real{0} && y<t_real{1} &&
		z>=t_real{0} && z<t_real{1};
}


/**
 * get triangle containing point pt
 */
template<class t_vec>
std::optional<std::size_t> get_containing_triag(const std::vector<std::vector<t_vec>>& triags,
	const t_vec& pt)
requires m::is_vec<t_vec>
{
	for(std::size_t idx=0; idx<triags.size(); ++idx)
	{
		const auto& triag = triags[idx];
		if(pt_inside_triag<t_vec>(triag[0], triag[1], triag[2], pt))
			return idx;
	}

	return std::nullopt;
}


/**
 * is delaunay triangle conflicting with point pt
 */
template<class t_vec>
bool is_conflicting_triag(const std::vector<t_vec>& triag, const t_vec& pt)
requires m::is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	// circumscribed circle radius
	t_vec center = calc_circumcentre<t_vec>(triag);
	t_real rad = m::norm<t_vec>(triag[0] - center);
	t_real dist = m::norm<t_vec>(pt - center);

	// point in circumscribed circle?
	return dist < rad;
}


/**
 * get delaunay triangles conflicting with point pt
 */
template<class t_vec>
std::vector<std::size_t> get_conflicting_triags(const std::vector<std::vector<t_vec>>& triags,
	const t_vec& pt)
requires m::is_vec<t_vec>
{
	std::vector<std::size_t> indices;

	for(std::size_t idx=0; idx<triags.size(); ++idx)
	{
		if(is_conflicting_triag<t_vec>(triags[idx], pt))
			indices.push_back(idx);
	}

	return indices;
}


template<class t_vec, class t_real = typename t_vec::value_type>
std::vector<t_vec>
_remove_duplicates(const std::vector<t_vec>& _verts, t_real eps=std::numeric_limits<t_real>::epsilon())
requires m::is_vec<t_vec>
{
	std::vector<t_vec> verts = _verts;

	// remove duplicate points
	verts.erase(std::unique(verts.begin(), verts.end(),
		[eps](const t_vec& vec1, const t_vec& vec2)->bool
		{ return m::equals<t_vec>(vec1, vec2, eps); }
		), verts.end());

	return verts;
}


template<class t_vec, class t_real = typename t_vec::value_type>
std::vector<t_vec>
_sort_vertices(const std::vector<t_vec>& _verts, t_real eps=std::numeric_limits<t_real>::epsilon())
requires m::is_vec<t_vec>
{
	std::vector<t_vec> verts = _verts;

	std::stable_sort(verts.begin(), verts.end(), [eps](const t_vec& vert1, const t_vec& vert2) -> bool
	{
		if(m::equals<t_real>(vert1[0], vert2[0], eps))
			return vert1[1] < vert2[1];
		return vert1[0] < vert2[0];
	});

	// remove unnecessary points
	auto iterCurX = verts.begin();
	for(auto iter = iterCurX; iter != verts.end(); std::advance(iter, 1))
	{
		// next x?
		if(!m::equals<t_real>((*iter)[0], (*iterCurX)[0], eps))
		{
			std::size_t num_same_x = iter - iterCurX;
			if(num_same_x >= 3)
				iter = verts.erase(std::next(iterCurX, 1), std::prev(iter, 1));
			iterCurX = iter;
		}
	}

	return verts;
}


template<class t_vec, class t_real = typename t_vec::value_type>
std::tuple<std::vector<t_vec>, t_vec>
sort_vertices_by_angle(const std::vector<t_vec>& _verts)
requires m::is_vec<t_vec>
{
	std::vector<t_vec> verts = _verts;
	//verts = _remove_duplicates<t_vec>(verts, eps);

	// sort by angle
	t_vec mean = std::accumulate(verts.begin(), verts.end(), m::zero<t_vec>(2));
	mean /= t_real(verts.size());
	std::stable_sort(verts.begin(), verts.end(), [&mean](const t_vec& vec1, const t_vec& vec2)->bool
	{ return line_angle<t_vec>(mean, vec1) < line_angle<t_vec>(mean, vec2); });

	return std::make_tuple(verts, mean);
}


/**
 * geometric series
 * see: https://en.wikipedia.org/wiki/Geometric_series
 */
template<class t_real = double>
t_real geo_series(t_real x, std::size_t n)
{
	/*t_real sum = 0.;
	 *	for(std::size_t i=0; i<n+1; ++i)
	 *		sum += std::pow(x, t_real(i));
	 *	std::cout << "sum: " << sum << std::endl;*/

	if(m::equals<t_real>(x, 1))
		return t_real(n+1) * x;
	else
		return (t_real(1) - std::pow(x, t_real(n+1))) / (t_real(1) - x);
}


template<class t_vec> requires m::is_vec<t_vec>
std::ostream& print_point(std::ostream& ostr, const t_vec& pt)
{
	ostr << "(" << pt[0] << ", " << pt[1] << ")";
	return ostr;
};


template<class t_vec, class t_line = std::pair<t_vec, t_vec>>
requires m::is_vec<t_vec>
std::ostream& print_line(std::ostream& ostr, const t_line& line)
{
	const auto& pt0 = std::get<0>(line);
	const auto& pt1 = std::get<1>(line);

	ostr << "(" << pt0[0] << ", " << pt0[1] << "), ("
	<< pt1[0] << ", " << pt1[1] << ")";

	return ostr;
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// convex hull algorithms
// @see (Klein 2005), ch. 4.1, pp. 155f
// @see (FUH 2020), ch. 3, pp. 113-160
// ----------------------------------------------------------------------------

/**
 * recursive calculation of convex hull
 * @see (FUH 2020), ch. 3.1.4, pp. 123-125
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::vector<t_vec>
_calc_hull_recursive_sorted(const std::vector<t_vec>& verts, t_real eps = 1e-5)
requires m::is_vec<t_vec>
{
	using namespace m_ops;

	// trivial cases to end recursion
	if(verts.size() <= 3)
	{
		std::vector<t_vec> hullverts;
		for(std::size_t vertidx=0; vertidx<verts.size(); ++vertidx)
			hullverts.push_back(verts[vertidx]);

		return std::get<0>(sort_vertices_by_angle<t_vec>(hullverts));
	}

	// divide
	std::size_t div = verts.size()/2;
	if(m::equals<t_real>(verts[div-1][0], verts[div][0], eps))
		++div;
	std::vector<t_vec> vertsLeft(verts.begin(), std::next(verts.begin(), div));
	std::vector<t_vec> vertsRight(std::next(verts.begin(), div), verts.end());

	// recurse
	std::vector<t_vec> hullLeft = _calc_hull_recursive_sorted(vertsLeft);
	std::vector<t_vec> hullRight = _calc_hull_recursive_sorted(vertsRight);


	// merge
	// upper part
	bool leftIsOnMax=false, rightIsOnMin=false;
	{
		auto _iterLeftMax = std::max_element(hullLeft.begin(), hullLeft.end(), [](const t_vec& vec1, const t_vec& vec2)->bool
		{ return vec1[0] < vec2[0]; });
		auto _iterRightMin = std::min_element(hullRight.begin(), hullRight.end(), [](const t_vec& vec1, const t_vec& vec2)->bool
		{ return vec1[0] < vec2[0]; });

		circular_wrapper circhullLeft(hullLeft);
		circular_wrapper circhullRight(hullRight);
		auto iterLeftMax = circhullLeft.begin() + (_iterLeftMax-hullLeft.begin());
		auto iterRightMin = circhullRight.begin() + (_iterRightMin-hullRight.begin());

		auto iterLeft = iterLeftMax;
		auto iterRight = iterRightMin;

		while(true)
		{
			bool leftChanged = false;
			bool rightChanged = false;

			while(side_of_line<t_vec>(*iterLeft, *iterRight, *(iterLeft+1)) > 0.)
			{
				++iterLeft;
				leftChanged = true;
			}
			while(side_of_line<t_vec>(*iterLeft, *iterRight, *(iterRight-1)) > 0.)
			{
				--iterRight;
				rightChanged = true;
			}

			// no more changes
			if(!leftChanged && !rightChanged)
				break;
		}

		if(iterLeft == iterLeftMax)
			leftIsOnMax = true;
		if(iterRight == iterRightMin)
			rightIsOnMin = true;

		circhullLeft.erase(iterLeftMax+1, iterLeft);
		circhullRight.erase(iterRight+1, iterRightMin);
	}

	// lower part
	{
		auto _iterLeftMax = std::max_element(hullLeft.begin(), hullLeft.end(), [](const t_vec& vec1, const t_vec& vec2)->bool
		{ return vec1[0] < vec2[0]; });
		auto _iterRightMin = std::min_element(hullRight.begin(), hullRight.end(), [](const t_vec& vec1, const t_vec& vec2)->bool
		{ return vec1[0] < vec2[0]; });

		circular_wrapper circhullLeft(hullLeft);
		circular_wrapper circhullRight(hullRight);
		auto iterLeftMax = circhullLeft.begin() + (_iterLeftMax-hullLeft.begin());
		auto iterRightMin = circhullRight.begin() + (_iterRightMin-hullRight.begin());

		auto iterLeft = iterLeftMax;
		auto iterRight = iterRightMin;

		while(true)
		{
			bool leftChanged = false;
			bool rightChanged = false;

			while(side_of_line<t_vec>(*iterLeft, *iterRight, *(iterLeft-1)) < 0.)
			{
				--iterLeft;
				leftChanged = true;
			}
			while(side_of_line<t_vec>(*iterLeft, *iterRight, *(iterRight+1)) < 0.)
			{
				++iterRight;
				rightChanged = true;
			}

			// no more changes
			if(!leftChanged && !rightChanged)
				break;
		}

		circhullLeft.erase(iterLeft+1, leftIsOnMax ? iterLeftMax : iterLeftMax+1);
		circhullRight.erase(rightIsOnMin ? iterRightMin+1 : iterRightMin, iterRight);
	}

	hullLeft.insert(hullLeft.end(), hullRight.begin(), hullRight.end());
	return std::get<0>(sort_vertices_by_angle<t_vec>(hullLeft));
}



template<class t_vec, class t_real = typename t_vec::value_type>
std::vector<t_vec>
calc_hull_recursive(const std::vector<t_vec>& _verts, t_real eps = 1e-5)
requires m::is_vec<t_vec>
{
	std::vector<t_vec> verts = _sort_vertices<t_vec>(_verts, eps);

	return _calc_hull_recursive_sorted<t_vec>(verts);
}


// ----------------------------------------------------------------------------


/**
 * tests if the vertex is in the hull
 */
template<class t_vec>
std::tuple<bool, std::size_t, std::size_t>
is_vert_in_hull(const std::vector<t_vec>& hull, const t_vec& newvert, const t_vec *vert_in_hull=nullptr)
requires m::is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	// get a point inside the hull if none given
	t_vec mean;
	if(!vert_in_hull)
	{
		mean = std::accumulate(hull.begin(), hull.end(), m::zero<t_vec>(2));
		mean /= t_real(hull.size());
		vert_in_hull = &mean;
	}

	for(std::size_t hullvertidx1=0; hullvertidx1<hull.size(); ++hullvertidx1)
	{
		std::size_t hullvertidx2 = hullvertidx1+1;
		if(hullvertidx2 >= hull.size())
			hullvertidx2 = 0;

		const t_vec& hullvert1 = hull[hullvertidx1];
		const t_vec& hullvert2 = hull[hullvertidx2];

		// new vertex is between these two points
		if(side_of_line<t_vec>(*vert_in_hull, hullvert1, newvert) > 0. &&
			side_of_line<t_vec>(*vert_in_hull, hullvert2, newvert) <= 0.)
		{
			// outside hull?
			if(side_of_line<t_vec>(hullvert1, hullvert2, newvert) < 0.)
				return std::make_tuple(false, hullvertidx1, hullvertidx2);
		}
	}
	return std::make_tuple(true, 0, 0);
};


/**
 * iterative calculation of convex hull
 * @see (FUH 2020), ch. 3.1.3, pp. 117-123
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::vector<t_vec>
calc_hull_iterative(const std::vector<t_vec>& _verts, t_real eps = 1e-5)
requires m::is_vec<t_vec>
{
	using namespace m_ops;

	std::vector<t_vec> verts = _remove_duplicates<t_vec>(_verts, eps);

	if(verts.size() <= 3)
		return verts;

	std::vector<t_vec> hull = {{ verts[0], verts[1], verts[2] }};
	t_vec vert_in_hull = m::zero<t_vec>(2);
	std::tie(hull, vert_in_hull) = sort_vertices_by_angle<t_vec>(hull);


	// insert new vertex into hull
	for(std::size_t vertidx=3; vertidx<verts.size(); ++vertidx)
	{
		const t_vec& newvert = verts[vertidx];

		// is the vertex already in the hull?
		auto [already_in_hull, hullvertidx1, hullvertidx2] =
			is_vert_in_hull<t_vec>(hull, newvert, &vert_in_hull);
		if(already_in_hull)
			continue;

		circular_wrapper circularverts(hull);
		auto iterLower = circularverts.begin() + hullvertidx1;
		auto iterUpper = circularverts.begin() + hullvertidx2;

		// correct cycles
		if(hullvertidx1 > hullvertidx2 && iterLower.GetRound()==iterUpper.GetRound())
			iterUpper.SetRound(iterLower.GetRound()+1);

		for(; iterLower.GetRound()>=-2; --iterLower)
		{
			if(side_of_line<t_vec>(*iterLower, newvert, *(iterLower-1)) >= 0.)
				break;
		}

		for(; iterUpper.GetRound()<=2; ++iterUpper)
		{
			if(side_of_line<t_vec>(*iterUpper, newvert, *(iterUpper+1)) <= 0.)
				break;
		}

		auto iter = iterUpper;
		if(iterLower+1 < iterUpper)
			iter = circularverts.erase(iterLower+1, iterUpper);
		hull.insert(iter.GetIter(), newvert);
	}

	return hull;
}


/**
 * iterative calculation of convex hull
 * @see (FUH 2020), ch. 3.1.3, pp. 117-123
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::vector<t_vec>
calc_hull_iterative_bintree(const std::vector<t_vec>& _verts, t_real eps = 1e-5)
requires m::is_vec<t_vec>
{
	using namespace m_ops;
	namespace intr = boost::intrusive;

	std::vector<t_vec> verts = _remove_duplicates<t_vec>(_verts, eps);

	if(verts.size() <= 3)
		return verts;

	std::vector<t_vec> starthull = {{ verts[0], verts[1], verts[2] }};
	t_vec vert_in_hull = std::accumulate(starthull.begin(), starthull.end(), m::zero<t_vec>(2));
	vert_in_hull /= t_real(starthull.size());


	using t_hook = intr::bs_set_member_hook<intr::link_mode<intr::normal_link>>;

	struct t_node
	{
		t_vec vert;
		t_real angle{};
		t_hook _h{};

		t_node(const t_vec& center, const t_vec& vert) : vert{vert}, angle{line_angle<t_vec>(center, vert)}
		{}

		bool operator<(const t_node& e2) const
		{
			return this->angle < e2.angle;
		}
	};

	using t_tree = intr::bstree<t_node, intr::member_hook<t_node, decltype(t_node::_h), &t_node::_h>>;
	t_tree hull;
	std::vector<t_node*> node_mem {{
		new t_node(vert_in_hull, verts[0]),
		new t_node(vert_in_hull, verts[1]),
		new t_node(vert_in_hull, verts[2]),
	}};

	for(t_node* node : node_mem)
		hull.insert_equal(*node);


	// test if the vertex is already in the hull
	auto is_in_hull = [&vert_in_hull, &hull](const t_vec& newvert)
		-> std::tuple<bool, std::size_t, std::size_t>
	{
		t_node tosearch(vert_in_hull, newvert);
		auto iter2 = hull.upper_bound(tosearch);
		// wrap around
		if(iter2 == hull.end())
			iter2 = hull.begin();

		auto iter1 = (iter2==hull.begin() ? std::next(hull.rbegin(),1).base() : std::prev(iter2,1));

		/*if(tosearch.angle < iter1->angle || tosearch.angle > iter2->angle)
			std::cerr << "angle: " << tosearch.angle/M_PI*180. << ", line1: " << iter1->angle/M_PI*180. << ", line2: " << iter2->angle/M_PI*180. << std::endl;*/

		const t_vec& vert1 = iter1->vert;
		const t_vec& vert2 = iter2->vert;

		// outside hull?
		if(side_of_line<t_vec>(vert1, vert2, newvert) < 0.)
		{
			std::size_t vertidx1 = std::distance(hull.begin(), iter1);
			std::size_t vertidx2 = std::distance(hull.begin(), iter2);

			return std::make_tuple(false, vertidx1, vertidx2);
		}

		return std::make_tuple(true, 0, 0);
	};


	// insert new vertex into hull
	for(std::size_t vertidx=3; vertidx<verts.size(); ++vertidx)
	{
		const t_vec& newvert = verts[vertidx];
		auto [already_in_hull, hullvertidx1, hullvertidx2] = is_in_hull(newvert);
		if(already_in_hull)
			continue;

		circular_wrapper circularverts(hull);
		auto iterLower = circularverts.begin()+hullvertidx1;
		auto iterUpper = circularverts.begin()+hullvertidx2;

		// correct cycles
		if(hullvertidx1 > hullvertidx2 && iterLower.GetRound()==iterUpper.GetRound())
			iterUpper.SetRound(iterLower.GetRound()+1);

		for(; iterLower.GetRound()>=-2; --iterLower)
		{
			if(side_of_line<t_vec>(iterLower->vert, newvert, (iterLower-1)->vert) >= 0.)
				break;
		}

		for(; iterUpper.GetRound()<=2; ++iterUpper)
		{
			if(side_of_line<t_vec>(iterUpper->vert, newvert, (iterUpper+1)->vert) <= 0.)
				break;
		}

		auto iter = iterUpper;
		if(std::distance(iterLower+1, iterUpper) > 0)
			iter = circularverts.erase(iterLower+1, iterUpper);

		t_node* newnode = new t_node(vert_in_hull, newvert);
		node_mem.push_back(newnode);
		hull.insert_equal(iter.GetIter(), *newnode);
	}


	// cleanups
	std::vector<t_vec> finalhull;
	for(auto iter = hull.begin(); iter != hull.end();)
	{
		finalhull.push_back(iter->vert);
		//t_node* node = &*iter;
		iter = hull.erase(iter);
		//delete node;
	}

	for(t_node* node : node_mem)
		delete node;

	return finalhull;
}


// ----------------------------------------------------------------------------


/**
 * calculation of convex hull
 * @see (FUH 2020), ch. 3.1.5, pp. 125-128
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::vector<t_vec>
calc_hull_contour(const std::vector<t_vec>& _verts, t_real eps = 1e-5)
requires m::is_vec<t_vec>
{
	using namespace m_ops;
	std::vector<t_vec> verts = _sort_vertices<t_vec>(_verts, eps);


	// contour determination
	{
		std::list<t_vec> contour_left_top, contour_left_bottom;
		std::pair<t_real, t_real> minmax_y_left
			= std::make_pair(std::numeric_limits<t_real>::max(), -std::numeric_limits<t_real>::max());

		for(const t_vec& vec : verts)
		{
			if(vec[1] > std::get<1>(minmax_y_left))
			{
				std::get<1>(minmax_y_left) = vec[1];
				contour_left_top.push_back(vec);
			}
			if(vec[1] < std::get<0>(minmax_y_left))
			{
				std::get<0>(minmax_y_left) = vec[1];
				contour_left_bottom.push_front(vec);
			}
		}


		std::list<t_vec> contour_right_top, contour_right_bottom;
		std::pair<t_real, t_real> minmax_y_right
			= std::make_pair(std::numeric_limits<t_real>::max(), -std::numeric_limits<t_real>::max());

		for(auto iter = verts.rbegin(); iter != verts.rend(); std::advance(iter, 1))
		{
			const t_vec& vec = *iter;
			if(vec[1] > std::get<1>(minmax_y_right))
			{
				std::get<1>(minmax_y_right) = vec[1];
				contour_right_top.push_front(vec);
			}
			if(vec[1] < std::get<0>(minmax_y_right))
			{
				std::get<0>(minmax_y_right) = vec[1];
				contour_right_bottom.push_back(vec);
			}
		}

		// convert to vector, only insert vertex if it's different than the last one
		verts.clear();
		verts.reserve(contour_left_top.size() + contour_right_top.size() +
			contour_left_bottom.size() + contour_right_bottom.size());
		for(const t_vec& vec : contour_left_top)
			if(!m::equals<t_vec>(*verts.rbegin(), vec, eps))
				verts.push_back(vec);
		for(const t_vec& vec : contour_right_top)
			if(!m::equals<t_vec>(*verts.rbegin(), vec, eps))
				verts.push_back(vec);
		for(const t_vec& vec : contour_right_bottom)
			if(!m::equals<t_vec>(*verts.rbegin(), vec, eps))
				verts.push_back(vec);
		for(const t_vec& vec : contour_left_bottom)
			if(!m::equals<t_vec>(*verts.rbegin(), vec, eps))
				verts.push_back(vec);

		if(verts.size() >= 2 && m::equals<t_vec>(*verts.begin(), *verts.rbegin(), eps))
			verts.erase(std::prev(verts.end(),1));
	}


	// hull calculation
	circular_wrapper circularverts(verts);
	for(std::size_t curidx = 1; curidx < verts.size()*2-1;)
	{
		if(curidx < 1)
			break;
		bool removed_points = false;

		// test convexity
		if(side_of_line<t_vec>(circularverts[curidx-1], circularverts[curidx+1], circularverts[curidx]) < 0.)
		{
			//std::cout << "vertex inside polygon: " << circularverts[curidx] << std::endl;
			for(std::size_t lastgood = curidx; lastgood >= 1; --lastgood)
			{
				if(side_of_line<t_vec>(circularverts[lastgood-1], circularverts[lastgood], circularverts[curidx+1]) <= 0.)
				{
					if(lastgood+1 > curidx+1)
						continue;

					circularverts.erase(std::next(circularverts.begin(), lastgood+1), std::next(circularverts.begin(), curidx+1));
					curidx = lastgood;
					removed_points = true;
					break;
				}
			}
		}

		if(!removed_points)
			++curidx;
	}

	return verts;
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// delaunay triangulation / voronoi regions
// @see (Klein 2005), ch. 6, pp. 269f
// @see (FUH 2020), ch. 5.3, pp. 228-232
// ----------------------------------------------------------------------------

/**
 * voronoi diagram for line segments
 * @see https://github.com/boostorg/polygon/blob/develop/example/voronoi_basic_tutorial.cpp
 * @see https://github.com/boostorg/polygon/blob/develop/example/voronoi_visual_utils.hpp
 * @see https://github.com/boostorg/polygon/blob/develop/example/voronoi_visualizer.cpp
 * @see https://www.boost.org/doc/libs/1_75_0/libs/polygon/doc/voronoi_diagram.htm
 */
template<class t_vec, class t_line=std::pair<t_vec, t_vec>, class t_graph=adjacency_matrix<typename t_vec::value_type>>
std::tuple<std::vector<t_vec>, std::vector<t_line>, std::vector<std::vector<t_vec>>, t_graph>
calc_voro(const std::vector<t_line>& lines)
requires m::is_vec<t_vec> && is_graph<t_graph>
{
#ifdef __GEO2D_USE_BOOST_POLY__
	using t_real = typename t_vec::value_type;
	namespace poly = boost::polygon;

	t_real parabola_eps = 1e-3;

	// length of infinite edges
	t_real infline_len = 1.;
	for(const t_line& line : lines)
	{
		t_vec dir = std::get<1>(line) - std::get<0>(line);
		t_real len = m::norm(dir);
		infline_len = std::max(infline_len, len);
	}
	infline_len *= 10.;

	using t_vorotraits = poly::voronoi_diagram_traits<t_real>;
	poly::voronoi_diagram<t_real, t_vorotraits> voro;
	poly::construct_voronoi(lines.begin(), lines.end(), &voro);

	// graph of voronoi vertices
	t_graph graph;

	//vertices
	std::vector<t_vec> vertices;
	for(const auto& vert : voro.vertices())
	{
		vertices.emplace_back(m::create<t_vec>({ vert.x(), vert.y() }));
		graph.AddVertex(std::to_string(vertices.size()));
	}


	auto get_vertex_idx = [&voro](const typename t_vorotraits::vertex_type* vert) -> std::optional<std::size_t>
	{
		// infinite edge?
		if(!vert)
			return std::nullopt;

		std::size_t idx = 0;
		for(const auto& vertex : voro.vertices())
		{
			if(&vertex == vert)
				return idx;
			++idx;
		}

		return std::nullopt;
	};


	// edges
	std::vector<std::vector<t_vec>> all_parabolic_edges;
	std::vector<t_line> linear_edges;
	linear_edges.reserve(voro.edges().size());

	for(const auto& edge : voro.edges())
	{
		// only bisectors, no internal edges
		if(edge.is_secondary())
			continue;

		// add graph edges
		const auto* vert0 = edge.vertex0();
		const auto* vert1 = edge.vertex1();
		auto vert0idx = get_vertex_idx(vert0);
		auto vert1idx = get_vertex_idx(vert1);
		if(vert0idx && vert1idx)
		{
			// TODO: arc length of parabolic edges
			t_real len = m::norm(vertices[*vert1idx] - vertices[*vert0idx]);
			graph.AddEdge(*vert0idx, *vert1idx, len);
			graph.AddEdge(*vert1idx, *vert0idx, len);
		}

		// get line segment
		auto get_segment = [&edge, &lines](bool twin) -> const t_line*
		{
			const auto* cell = twin ? edge.twin()->cell() : edge.cell();
			if(!cell)
				return nullptr;

			const t_line& line = lines[cell->source_index()];
			return &line;
		};


		// get line segment endpoint
		auto get_segment_point = [&edge, &get_segment](bool twin) -> const t_vec*
		{
			const auto* cell = twin ? edge.twin()->cell() : edge.cell();
			if(!cell)
				return nullptr;

			const t_line* line = get_segment(twin);
			const t_vec* vec = nullptr;
			if(!line)
				return nullptr;

			switch(cell->source_category())
			{
				case poly::SOURCE_CATEGORY_SEGMENT_START_POINT:
					vec = &std::get<0>(*line);
					break;
				case poly::SOURCE_CATEGORY_SEGMENT_END_POINT:
					vec = &std::get<1>(*line);
					break;
				default:
					break;
			}

			return vec;
		};


		// converter functions
		auto to_point_data = [](const t_vec& vec) -> poly::point_data<t_real>
		{
			return poly::point_data<t_real>{vec[0], vec[1]};
		};

		auto vertex_to_point_data = [](const typename t_vorotraits::vertex_type& vec) -> poly::point_data<t_real>
		{
			return poly::point_data<t_real>{vec.x(), vec.y()};
		};

		auto to_vec = [](const poly::point_data<t_real>& pt) -> t_vec
		{
			return m::create<t_vec>({ pt.x(), pt.y() });
		};

		auto vertex_to_vec = [](const typename t_vorotraits::vertex_type& vec) -> t_vec
		{
			return m::create<t_vec>({ vec.x(), vec.y() });
		};

		auto to_segment_data = [&to_point_data](const t_line& line) -> poly::segment_data<t_real>
		{
			auto pt1 = to_point_data(std::get<0>(line));
			auto pt2 = to_point_data(std::get<1>(line));

			return poly::segment_data<t_real>{pt1, pt2};
		};


		// parabolic edge
		if(edge.is_curved())
		{
			const t_line* seg = get_segment(edge.cell()->contains_point());
			const t_vec* pt = get_segment_point(!edge.cell()->contains_point());
			if(!seg || !pt)
				continue;

			std::vector<poly::point_data<t_real>> parabola
				{{ vertex_to_point_data(*edge.vertex0()), vertex_to_point_data(*edge.vertex1()) }};

			poly::voronoi_visual_utils<t_real>::discretize(
				to_point_data(*pt), to_segment_data(*seg),
				parabola_eps, &parabola);

			std::vector<t_vec> parabolic_edges;
			parabolic_edges.reserve(parabola.size());
			for(const auto& parabola_pt : parabola)
				parabolic_edges.emplace_back(to_vec(parabola_pt));
			all_parabolic_edges.emplace_back(std::move(parabolic_edges));
		}

		// linear edge
		else
		{
			// finite edge
			if(edge.is_finite())
			{
				linear_edges.push_back(std::make_pair(
					vertex_to_vec(*edge.vertex0()), vertex_to_vec(*edge.vertex1())));
			}

			// infinite edge
			else
			{
				t_vec lineorg;
				bool inverted = false;
				if(edge.vertex0())
				{
					lineorg = vertex_to_vec(*edge.vertex0());
					inverted = false;
				}
				else if(edge.vertex1())
				{
					lineorg = vertex_to_vec(*edge.vertex1());
					inverted = true;
				}
				else
				{
					continue;
				}

				const t_vec* vec = get_segment_point(false);
				const t_vec* twinvec = get_segment_point(true);

				if(!vec || !twinvec)
					continue;

				t_vec perpdir = *vec - *twinvec;
				if(inverted)
					perpdir = -perpdir;
				t_vec linedir = m::create<t_vec>({ perpdir[1], -perpdir[0] });

				linedir /= m::norm(linedir);
				linedir *= infline_len;

				linear_edges.push_back(std::make_pair(lineorg, lineorg + linedir));
			}
		}
	}

	return std::make_tuple(vertices, linear_edges, all_parabolic_edges, graph);

#else

	// disable function
	return std::make_tuple(std::vector<t_vec>{},
		std::vector<t_line>{}, std::vector<std::vector<t_vec>>{},
		t_graph{});

#endif
}



/**
 * delaunay triangulation and voronoi vertices
 * @returns [ voronoi vertices, triangles, neighbour triangle indices ]
 */
template<class t_vec>
std::tuple<std::vector<t_vec>, std::vector<std::vector<t_vec>>, std::vector<std::set<std::size_t>>>
calc_delaunay(int dim, const std::vector<t_vec>& verts, bool only_hull)
requires m::is_vec<t_vec>
{
	//using namespace m_ops;
	namespace qh = orgQhull;

	using t_real = typename t_vec::value_type;
	using t_real_qhull = coordT;

	std::vector<t_vec> voronoi;						// voronoi vertices
	std::vector<std::vector<t_vec>> triags;			// delaunay triangles
	std::vector<std::set<std::size_t>> neighbours;	// neighbour triangle indices

	try
	{
		std::vector<t_real_qhull> _verts;
		_verts.reserve(verts.size() * dim);
		for(const t_vec& vert : verts)
			for(int i=0; i<dim; ++i)
				_verts.push_back(t_real_qhull{vert[i]});

		qh::Qhull qh{"triag", dim, int(_verts.size()/dim), _verts.data(), only_hull ? "Qt" : "v Qu QJ" };
		if(qh.hasQhullMessage())
			std::cout << qh.qhullMessage() << std::endl;


		//qh::QhullVertexList vertices{qh.vertexList()};
		qh::QhullFacetList facets{qh.facetList()};
		std::vector<void*> facetHandles{};

		// get all triangles
		for(auto iterFacet=facets.begin(); iterFacet!=facets.end(); ++iterFacet)
		{
			if(iterFacet->isUpperDelaunay())
				continue;
			facetHandles.push_back(iterFacet->getBaseT());

			if(!only_hull)
			{
				qh::QhullPoint pt = iterFacet->voronoiVertex();

				t_vec vec = m::create<t_vec>(dim);
				for(int i=0; i<dim; ++i)
					vec[i] = t_real{pt[i]};

				voronoi.emplace_back(std::move(vec));
			}


			std::vector<t_vec> thetriag;
			qh::QhullVertexSet vertices = iterFacet->vertices();

			for(auto iterVertex=vertices.begin(); iterVertex!=vertices.end(); ++iterVertex)
			{
				qh::QhullPoint pt = (*iterVertex).point();

				t_vec vec = m::create<t_vec>(dim);
				for(int i=0; i<dim; ++i)
					vec[i] = t_real{pt[i]};

				thetriag.emplace_back(std::move(vec));
			}

			if(dim == 2)
				std::tie(thetriag, std::ignore) = sort_vertices_by_angle<t_vec>(thetriag);
			triags.emplace_back(std::move(thetriag));
		}


		// find neighbouring triangles
		if(!only_hull)
		{
			neighbours.resize(triags.size());

			std::size_t facetIdx = 0;
			for(auto iterFacet=facets.begin(); iterFacet!=facets.end(); ++iterFacet)
			{
				if(iterFacet->isUpperDelaunay())
					continue;

				qh::QhullFacetSet neighbourFacets{iterFacet->neighborFacets()};
				for(auto iterNeighbour=neighbourFacets.begin(); iterNeighbour!=neighbourFacets.end(); ++iterNeighbour)
				{
					void* handle = (*iterNeighbour).getBaseT();
					auto iterHandle = std::find(facetHandles.begin(), facetHandles.end(), handle);
					if(iterHandle != facetHandles.end())
					{
						std::size_t handleIdx = iterHandle - facetHandles.begin();
						neighbours[facetIdx].insert(handleIdx);
					}
				}

				if(++facetIdx >= triags.size())
					break;
			}
		}
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}

	return std::make_tuple(voronoi, triags, neighbours);
}



/**
 * @returns [triangle index, shared index 1, shared index 2, non-shared index]
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::optional<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>
get_triag_sharing_edge(std::vector<std::vector<t_vec>>& triags,
	const t_vec& vert1, const t_vec& vert2, std::size_t curtriagidx, t_real eps = 1e-5)
requires m::is_vec<t_vec>
{
	for(std::size_t i=0; i<triags.size(); ++i)
	{
		if(i == curtriagidx)
			continue;

		const auto& triag = triags[i];

		// test all edge combinations
		if(m::equals<t_vec>(triag[0], vert1, eps) && m::equals<t_vec>(triag[1], vert2, eps))
			return std::make_tuple(i, 0, 1, 2);
		if(m::equals<t_vec>(triag[1], vert1, eps) && m::equals<t_vec>(triag[0], vert2, eps))
			return std::make_tuple(i, 1, 0, 2);
		if(m::equals<t_vec>(triag[0], vert1, eps) && m::equals<t_vec>(triag[2], vert2, eps))
			return std::make_tuple(i, 0, 2, 1);
		if(m::equals<t_vec>(triag[2], vert1, eps) && m::equals<t_vec>(triag[0], vert2, eps))
			return std::make_tuple(i, 2, 0, 1);
		if(m::equals<t_vec>(triag[1], vert1, eps) && m::equals<t_vec>(triag[2], vert2, eps))
			return std::make_tuple(i, 1, 2, 0);
		if(m::equals<t_vec>(triag[2], vert1, eps) && m::equals<t_vec>(triag[1], vert2, eps))
			return std::make_tuple(i, 2, 1, 0);
	}

	// no shared edge found
	return std::nullopt;
}



template<class t_vec, class t_real = typename t_vec::value_type>
void flip_edge(std::vector<std::vector<t_vec>>& triags,
	std::size_t triagidx, std::size_t nonsharedidx, t_real eps = 1e-5)
requires m::is_vec<t_vec>
{
	std::size_t sharedidx1 = (nonsharedidx+1) % triags[triagidx].size();
	std::size_t sharedidx2 = (nonsharedidx+2) % triags[triagidx].size();

	// get triangle on other side of shared edge
	auto optother = get_triag_sharing_edge(
		triags, triags[triagidx][sharedidx1], triags[triagidx][sharedidx2], triagidx, eps);
	if(!optother)
		return;
	const auto [othertriagidx, othersharedidx1, othersharedidx2, othernonsharedidx] = *optother;

	if(is_conflicting_triag<t_vec>(triags[othertriagidx], triags[triagidx][nonsharedidx]))
	{
		triags[triagidx] = std::vector<t_vec>
		{{
			triags[triagidx][nonsharedidx],
			triags[othertriagidx][othernonsharedidx],
			triags[othertriagidx][othersharedidx1]
		}};

		triags[othertriagidx] = std::vector<t_vec>
		{{
			triags[triagidx][nonsharedidx],
			triags[othertriagidx][othernonsharedidx],
			triags[othertriagidx][othersharedidx2]
		}};

		// also check neighbours of newly created triangles for conflicts
		flip_edge(triags, othertriagidx, othernonsharedidx, eps);
		flip_edge(triags, othertriagidx, othersharedidx1, eps);
		flip_edge(triags, othertriagidx, othersharedidx2, eps);

		flip_edge(triags, triagidx, nonsharedidx, eps);
		flip_edge(triags, triagidx, sharedidx1, eps);
		flip_edge(triags, triagidx, sharedidx2, eps);
	}
}



/**
 * iterative delaunay triangulation
 * @see (FUH 2020), ch. 6.2, pp. 269-282
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::tuple<std::vector<t_vec>, std::vector<std::vector<t_vec>>, std::vector<std::set<std::size_t>>>
calc_delaunay_iterative(const std::vector<t_vec>& verts , t_real eps = 1e-5)
requires m::is_vec<t_vec>
{
	using namespace m_ops;

	std::vector<t_vec> voronoi;						// voronoi vertices
	std::vector<std::vector<t_vec>> triags;			// delaunay triangles
	std::vector<std::set<std::size_t>> neighbours;	// neighbour triangle indices

	if(verts.size() < 3)
		return std::make_tuple(voronoi, triags, neighbours);

	// first triangle
	triags.emplace_back(std::vector<t_vec>{{ verts[0], verts[1], verts[2] }});

	// currently inserted vertices
	std::vector<t_vec> curverts{{ verts[0], verts[1], verts[2] }};

	// insert vertices iteratively
	for(std::size_t newvertidx=3; newvertidx<verts.size(); ++newvertidx)
	{
		const t_vec& newvert = verts[newvertidx];
		//std::cout << "newvert " << newvertidx-3 << ": " << newvert << std::endl;

		// find triangle containing the new vertex
		if(auto optidx = get_containing_triag<t_vec>(triags, newvert); optidx)
		{
			//std::cout << "inside" << std::endl;

			auto conttriag = std::move(triags[*optidx]);
			triags.erase(triags.begin() + *optidx);

			// new delaunay edges connecting to newvert
			triags.emplace_back(std::vector<t_vec>{{ newvert, conttriag[0], conttriag[1] }});
			triags.emplace_back(std::vector<t_vec>{{ newvert, conttriag[0], conttriag[2] }});
			triags.emplace_back(std::vector<t_vec>{{ newvert, conttriag[1], conttriag[2] }});

			flip_edge(triags, triags.size()-3, 0, eps);
			flip_edge(triags, triags.size()-2, 0, eps);
			flip_edge(triags, triags.size()-1, 0, eps);
		}

		// new vertex is outside of any triangle
		else
		{
			//std::cout << "outside" << std::endl;
			auto hull = calc_hull_iterative_bintree<t_vec>(curverts, eps);
			std::tie(hull, std::ignore) = sort_vertices_by_angle<t_vec>(hull);

			// find the points in the hull visible from newvert
			std::vector<t_vec> visible;
			{
				// start indices
				auto [already_in_hull, hullvertidx1, hullvertidx2] =
					is_vert_in_hull<t_vec>(hull, newvert);
				if(already_in_hull)
					continue;

				// find visible vertices like in calc_hull_iterative
				circular_wrapper circularverts(hull);
				auto iterLower = circularverts.begin() + hullvertidx1;
				auto iterUpper = circularverts.begin() + hullvertidx2;

				// correct cycles
				if(hullvertidx1 > hullvertidx2 && iterLower.GetRound()==iterUpper.GetRound())
					iterUpper.SetRound(iterLower.GetRound()+1);

				for(; iterLower.GetRound()>=-2; --iterLower)
				{
					if(side_of_line<t_vec>(*iterLower, newvert, *(iterLower-1)) >= 0.)
						break;
				}

				for(; iterUpper.GetRound()<=2; ++iterUpper)
				{
					if(side_of_line<t_vec>(*iterUpper, newvert, *(iterUpper+1)) <= 0.)
						break;
				}

				for(auto iter=iterLower; iter<=iterUpper; ++iter)
				{
					//std::cout << "vis: " << *iter << std::endl;
					visible.push_back(*iter);
				}
			}

			for(std::size_t visidx=0; visidx<visible.size()-1; ++visidx)
			{
				triags.emplace_back(std::vector<t_vec>{{ newvert, visible[visidx], visible[visidx+1] }});
				flip_edge(triags, triags.size()-1, 0, eps);
			}
		}

		//std::cout << "----------------------------------------" << std::endl;
		curverts.push_back(newvert);
	}


	// find neighbouring triangles and voronoi vertices
	neighbours.resize(triags.size());

	for(std::size_t triagidx=0; triagidx<triags.size(); ++triagidx)
	{
		// sort vertices
		auto& triag = triags[triagidx];
		std::tie(triag, std::ignore) = sort_vertices_by_angle<t_vec>(triag);

		// voronoi vertices
		voronoi.emplace_back(calc_circumcentre<t_vec>(triag));

		// neighbouring triangle indices
		auto optother1 = get_triag_sharing_edge(triags, triags[triagidx][0], triags[triagidx][1], triagidx, eps);
		auto optother2 = get_triag_sharing_edge(triags, triags[triagidx][0], triags[triagidx][2], triagidx, eps);
		auto optother3 = get_triag_sharing_edge(triags, triags[triagidx][1], triags[triagidx][2], triagidx, eps);

		if(optother1) neighbours[triagidx].insert(std::get<0>(*optother1));
		if(optother2) neighbours[triagidx].insert(std::get<0>(*optother2));
		if(optother3) neighbours[triagidx].insert(std::get<0>(*optother3));
	}


	return std::make_tuple(voronoi, triags, neighbours);
}



/**
 * delaunay triangulation using parabolic trafo
 * @see (Berg 2008), pp. 254-256 and p. 168
 * @see (FUH 2020), ch. 6.5, pp. 298-300
 */
template<class t_vec>
std::tuple<std::vector<t_vec>, std::vector<std::vector<t_vec>>, std::vector<std::set<std::size_t>>>
calc_delaunay_parabolic(const std::vector<t_vec>& verts)
requires m::is_vec<t_vec>
{
	using namespace m_ops;
	namespace qh = orgQhull;

	using t_real = typename t_vec::value_type;
	using t_real_qhull = coordT;

	const int dim = 2;
	std::vector<t_vec> voronoi;						// voronoi vertices
	std::vector<std::vector<t_vec>> triags;			// delaunay triangles
	std::vector<std::set<std::size_t>> neighbours;	// neighbour triangle indices

	try
	{
		std::vector<t_real_qhull> _verts;
		_verts.reserve(verts.size()*(dim+1));
		for(const t_vec& vert : verts)
		{
			_verts.push_back(t_real_qhull{vert[0]});
			_verts.push_back(t_real_qhull{vert[1]});
			_verts.push_back(t_real_qhull{vert[0]*vert[0] + vert[1]*vert[1]});
		}

		qh::Qhull qh{"triag", dim+1, int(_verts.size()/(dim+1)), _verts.data(), "Qt"};
		if(qh.hasQhullMessage())
			std::cout << qh.qhullMessage() << std::endl;


		qh::QhullFacetList facets{qh.facetList()};
		std::vector<void*> facetHandles{};


		auto facetAllowed = [](auto iterFacet) -> bool
		{
			if(iterFacet->isUpperDelaunay())
				return false;

			// filter out non-visible part of hull
			qh::QhullHyperplane plane = iterFacet->hyperplane();
			t_vec normal = m::create<t_vec>(dim+1);
			for(int i=0; i<dim+1; ++i)
				normal[i] = t_real{plane[i]};
			// normal pointing upwards?
			if(normal[2] > 0.)
				return false;

			return true;
		};


		for(auto iterFacet=facets.begin(); iterFacet!=facets.end(); ++iterFacet)
		{
			if(!facetAllowed(iterFacet))
				continue;

			std::vector<t_vec> thetriag;
			qh::QhullVertexSet vertices = iterFacet->vertices();

			for(auto iterVertex=vertices.begin(); iterVertex!=vertices.end(); ++iterVertex)
			{
				qh::QhullPoint pt = (*iterVertex).point();

				t_vec vec = m::create<t_vec>(dim);
				for(int i=0; i<dim; ++i)
					vec[i] = t_real{pt[i]};

				thetriag.emplace_back(std::move(vec));
			}

			voronoi.emplace_back(calc_circumcentre<t_vec>(thetriag));
			std::tie(thetriag, std::ignore) = sort_vertices_by_angle<t_vec>(thetriag);
			triags.emplace_back(std::move(thetriag));
			facetHandles.push_back(iterFacet->getBaseT());
		}


		// find neighbouring triangles
		neighbours.resize(triags.size());

		std::size_t facetIdx = 0;
		for(auto iterFacet=facets.begin(); iterFacet!=facets.end(); ++iterFacet)
		{
			if(!facetAllowed(iterFacet))
				continue;

			qh::QhullFacetSet neighbourFacets{iterFacet->neighborFacets()};
			for(auto iterNeighbour=neighbourFacets.begin(); iterNeighbour!=neighbourFacets.end(); ++iterNeighbour)
			{
				void* handle = (*iterNeighbour).getBaseT();
				auto iterHandle = std::find(facetHandles.begin(), facetHandles.end(), handle);
				if(iterHandle != facetHandles.end())
				{
					std::size_t handleIdx = iterHandle - facetHandles.begin();
					neighbours[facetIdx].insert(handleIdx);
				}
			}

			if(++facetIdx >= triags.size())
				break;
		}
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}

	return std::make_tuple(voronoi, triags, neighbours);
}



/**
 * get all edges from a delaunay triangulation
 */
template<class t_vec,
	class t_edge = std::pair<std::size_t, std::size_t>,
	class t_real = typename t_vec::value_type>
std::vector<t_edge>
get_edges(const std::vector<t_vec>& verts, const std::vector<std::vector<t_vec>>& triags, t_real eps)
{
	auto get_vert_idx = [&verts, eps](const t_vec& vert) -> std::optional<std::size_t>
	{
		for(std::size_t vertidx=0; vertidx<verts.size(); ++vertidx)
		{
			const t_vec& vert2 = verts[vertidx];
			if(m::equals<t_vec>(vert, vert2, eps))
				return vertidx;
		}

		return std::nullopt;
	};


	std::vector<t_edge> edges;

	for(std::size_t vertidx=0; vertidx<verts.size(); ++vertidx)
	{
		const t_vec& vert = verts[vertidx];

		for(const auto& triag : triags)
		{
			for(std::size_t i=0; i<triag.size(); ++i)
			{
				const t_vec& triagvert = triag[i];

				if(m::equals<t_vec>(vert, triagvert, eps))
				{
					const t_vec& vert2 = triag[(i+1) % triag.size()];
					const t_vec& vert3 = triag[(i+2) % triag.size()];

					std::size_t vert2idx = *get_vert_idx(vert2);
					std::size_t vert3idx = *get_vert_idx(vert3);

					edges.push_back(std::make_pair(vertidx, vert2idx));
					edges.push_back(std::make_pair(vertidx, vert3idx));

					//std::cout << vert2idx << " -> " << vertidx << " -> " << vert3idx << std::endl;
				}
			}
		}
	}

	return edges;
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// spanning tree
// ----------------------------------------------------------------------------

/**
 * finds loops in an undirected graph
 */
template<class t_edge = std::pair<std::size_t, std::size_t>>
bool has_loops(const std::vector<t_edge>& edges, std::size_t start_from, std::size_t start_to)
{
	// [from, to]
	std::stack<t_edge> tovisit;
	tovisit.push(std::make_pair(start_from, start_to));

	std::set<std::size_t> visitedverts;
	visitedverts.insert(start_from);

	std::set<t_edge> visitededges;

	// visit connected vertices
	while(!tovisit.empty())
	{
		auto topedge = tovisit.top();
		auto [vertfrom, vertto] = topedge;
		tovisit.pop();

		if(visitededges.find(topedge) != visitededges.end())
			continue;

		visitededges.insert(std::make_pair(vertfrom, vertto));
		visitededges.insert(std::make_pair(vertto, vertfrom));

		// has this vertex already been visited? => loop in graph
		if(visitedverts.find(vertto) != visitedverts.end())
			return true;

		visitedverts.insert(vertto);

		// get all edges from current vertex
		for(auto iter=edges.begin(); iter!=edges.end(); ++iter)
		{
			// forward direction
			if(iter->first == vertto)
				tovisit.push(std::make_pair(iter->first, iter->second));
			// backward direction
			if(iter->second == vertto)
				tovisit.push(std::make_pair(iter->second, iter->first));
		}
	}

	return false;
}


/**
 * minimal spanning tree
 * @see (FUH 2020), ch. 5.2.3, pp. 221-224
 * @see https://de.wikipedia.org/wiki/Algorithmus_von_Kruskal
 */
template<class t_vec, class t_edge = std::pair<std::size_t, std::size_t>>
std::vector<t_edge>
calc_min_spantree(const std::vector<t_vec>& verts, const std::vector<t_edge>& _edges)
requires m::is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	std::vector<t_edge> edges = _edges;

	std::stable_sort(edges.begin(), edges.end(), [&verts](const t_edge& edge1, const t_edge& edge2) -> bool
	{
		t_vec dir1 = verts[edge1.first]-verts[edge1.second];
		t_vec dir2 = verts[edge2.first]-verts[edge2.second];

		t_real len1sq = m::inner(dir1, dir1);
		t_real len2sq = m::inner(dir2, dir2);

		return len1sq >= len2sq;
	});


	std::vector<t_edge> span;

	while(edges.size())
	{
		t_edge edge = std::move(edges.back());
		edges.pop_back();

		span.push_back(edge);
		if(has_loops<t_edge>(span, edge.first, edge.second))
			span.pop_back();
	}

	return span;
}


/**
 * minimum spanning tree
 */
template<class t_vec, class t_edge = std::pair<std::size_t, std::size_t>>
std::vector<t_edge>
calc_min_spantree_boost(const std::vector<t_vec>& verts)
requires m::is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	struct t_edge_weight
	{
		t_edge_weight() = default;
		t_edge_weight(t_real weight) : weight(weight) {}

		t_real weight{1};
	};

	using t_graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, t_vec, t_edge_weight>;
	using t_edge_descr = typename boost::graph_traits<t_graph>::edge_descriptor;

	t_graph graph;
	auto weight = boost::get(&t_edge_weight::weight, graph);

	for(const t_vec& vert : verts)
		boost::add_vertex(vert, graph);

	for(std::size_t i=0; i<verts.size(); ++i)
	{
		const t_vec& vert1 = verts[i];
		for(std::size_t j=i+1; j<verts.size(); ++j)
		{
			const t_vec& vert2 = verts[j];
			t_real dist = m::norm(vert2-vert1);
			boost::add_edge(boost::vertex(i, graph), boost::vertex(j, graph), t_edge_weight{dist}, graph);
		}
	}

	std::vector<t_edge_descr> spanning_edges;
	boost::kruskal_minimum_spanning_tree(graph, std::back_inserter(spanning_edges), boost::weight_map(weight));

	std::vector<t_edge> span;
	for(auto iter=spanning_edges.begin(); iter!=spanning_edges.end(); ++iter)
	{
		std::size_t idx1 = boost::source(*iter, graph);
		std::size_t idx2 = boost::target(*iter, graph);
		span.emplace_back(std::make_pair(idx1, idx2));
	}

	return span;
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// kernel
// @see (Klein 2005), ch. 4.4, pp. 195f
// @see (FUH 2020), ch. 3.3.2, pp. 142-143
// ----------------------------------------------------------------------------

template<class t_vec, class t_real = typename t_vec::value_type>
std::pair<bool, t_vec> get_inters_dual(const t_vec& vert1, const t_vec& vert2)
requires m::is_vec<t_vec>
{
	t_real slope1 = -vert1[0];
	t_real offs1 = vert1[1];
	t_real slope2 = -vert2[0];
	t_real offs2 = vert2[1];

	// slope1*x + offs1 = slope2*x + offs2
	// x = (offs2 - offs1) / (slope1 - slope2)
	t_real intersX = (offs2-offs1) / (slope1-slope2);
	t_real intersY = offs1 + slope1*intersX;

	if(std::isfinite(intersX) && std::isfinite(intersY))
		return std::make_pair(true, m::create<t_vec>({intersX, intersY}));
	return std::make_pair(false, t_vec{});
};



/**
 * lower halfplane intersection vertices
 * @see (FUH 2020), ch. 3.3.2, pp. 142-143
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::vector<t_vec>
halfplaneverts(const std::vector<t_vec>& verts, t_real eps)
requires m::is_vec<t_vec>
{
	//using namespace m_ops;
	if(verts.size() < 3)
		return std::vector<t_vec>({});


	std::vector<t_vec> vertsDual;
	for(std::size_t vertidx=0; vertidx<verts.size(); ++vertidx)
	{
		std::size_t vertidxNext = (vertidx+1) % verts.size();

		const t_vec& vert1 = verts[vertidx];
		const t_vec& vert2 = verts[vertidxNext];

		if(t_vec dir=vert2-vert1; -dir[0] > 0.)
		{
			t_real slope = (vert2[1]-vert1[1]) / (vert2[0]-vert1[0]);
			t_real offs = vert1[1] - vert1[0]*slope;

			vertsDual.emplace_back(m::create<t_vec>({ -slope, offs }));
		}
	}


	std::vector<t_vec> hullDual = calc_hull_iterative_bintree<t_vec>(vertsDual, eps);
	std::tie(hullDual, std::ignore) = sort_vertices_by_angle<t_vec>(hullDual);


	std::vector<t_vec> intersverts;

	for(std::size_t hullidx=0; hullidx<hullDual.size(); ++hullidx)
	{
		std::size_t hullidx2 = (hullidx + 1) % hullDual.size();

		const t_vec& vec1 = hullDual[hullidx];
		const t_vec& vec2 = hullDual[hullidx2];

		t_vec dir = vec2 - vec1;
		//std::cout << "edge: " << vec1 << "\t" << vec2 << "\tdir: " << dir << std::endl;

		if(-dir[0] < 0.)
		{
			if(auto [ok, inters] = get_inters_dual(vec1, vec2); ok)
			{
				//std::cout << "intersection: " << inters << std::endl;
				intersverts.emplace_back(std::move(inters));
			}
		}
	}
	//std::cout << std::endl;

	std::tie(intersverts, std::ignore) = sort_vertices_by_angle<t_vec>(intersverts);
	return intersverts;
}


/**
 * kernel
 * @see (FUH 2020), ch. 3.3.2, pp. 142-143
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::vector<t_vec>
calc_ker_tst(const std::vector<t_vec>& verts, t_real eps)
requires m::is_vec<t_vec>
{
	if(verts.size() < 3)
		return std::vector<t_vec>({});

	std::vector<t_vec> vertsneg;
	for(const t_vec& vec : verts)
		vertsneg.emplace_back(-vec);

	auto verts1 = halfplaneverts(verts, eps);
	auto verts2 = halfplaneverts(vertsneg, eps);

	for(const t_vec& vec : verts2)
		verts1.emplace_back(-vec);

	std::tie(verts1, std::ignore) = sort_vertices_by_angle<t_vec>(verts1);
 	return verts1;
}



template<class t_vec, class t_edge = std::pair<t_vec, t_vec>, class t_real = typename t_vec::value_type>
std::vector<t_vec>
ker_from_edges(const std::vector<t_edge>& edges, t_real eps)
requires m::is_vec<t_vec>
{
	std::vector<t_vec> intersections;

	for(std::size_t i=0; i<edges.size(); ++i)
	{
		for(std::size_t j=i+1; j<edges.size(); ++j)
		{
			if(auto [ok, inters] = intersect_lines<t_vec>(
				edges[i].first, edges[i].second,
				edges[j].first, edges[j].second, false, eps); ok)
			{
				intersections.emplace_back(std::move(inters));
			}
		}
	}


	std::vector<t_vec> ker;

	for(const t_vec& inters : intersections)
	{
		bool in_ker = true;
		for(const t_edge& edge : edges)
		{
			if(side_of_line<t_vec>(edge.first, edge.second, inters) < -eps)
			{
				in_ker = false;
				break;
			}
		}

		if(in_ker)
			ker.push_back(inters);
	}

	std::tie(ker, std::ignore) = sort_vertices_by_angle<t_vec>(ker);
	return ker;
}



/**
 * kernel of a polygon
 * (still O(n^2)! TODO: only check contributing edges)
 * vertices have to be sorted in ccw order
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::vector<t_vec>
calc_ker(const std::vector<t_vec>& verts, t_real eps)
requires m::is_vec<t_vec>
{
	if(verts.size() < 3)
		return std::vector<t_vec>({});

	using t_edge = std::pair<t_vec, t_vec>;
	std::vector<t_edge> edges;

	for(std::size_t vertidx=0; vertidx<verts.size(); ++vertidx)
	{
		std::size_t vertidxNext = (vertidx+1) % verts.size();
		edges.emplace_back(std::make_pair(verts[vertidx], verts[vertidxNext]));
	}

	return ker_from_edges<t_vec, t_edge>(edges, eps);
}



// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// line segment intersections
// @see (Klein 2005), ch. 2.3.2, pp. 64f
// @see (FUH 2020), ch. 2.3.2, pp. 69-80
// ----------------------------------------------------------------------------

template<class t_vec, class t_line = std::pair<t_vec, t_vec>, class t_real = typename t_vec::value_type>
std::vector<std::tuple<std::size_t, std::size_t, t_vec>>
intersect_ineff(const std::vector<t_line>& lines, t_real eps = 1e-6)
requires m::is_vec<t_vec>
{
	std::vector<std::tuple<std::size_t, std::size_t, t_vec>> intersections;

	for(std::size_t i=0; i<lines.size(); ++i)
	{
		for(std::size_t j=i+1; j<lines.size(); ++j)
		{
			const t_line& line1 = lines[i];
			const t_line& line2 = lines[j];

			if(auto [intersects, pt] = intersect_lines<t_line>(line1, line2, eps); intersects)
				intersections.emplace_back(std::make_tuple(i, j, pt));
		}
	}

	return intersections;
}


template<class t_vec, class t_line = std::pair<t_vec, t_vec>>
requires m::is_vec<t_vec>
bool cmp_line(const t_line& line1, const t_line& line2,
	typename t_vec::value_type x, typename t_vec::value_type eps)
{
	using t_real = typename t_vec::value_type;

	//t_real line1_y = get_line_y<t_vec>(line1, x1);
	//t_real line2_y = get_line_y<t_vec>(line2, x2);

	auto [slope1, offs1] = get_line_slope_offs<t_vec>(line1);
	auto [slope2, offs2] = get_line_slope_offs<t_vec>(line2);
	t_real line1_y = slope1*x + offs1;
	t_real line2_y = slope2*x + offs2;

	// equal y -> compare by slope
	if(m::equals<t_real>(line1_y, line2_y, eps))
		return slope1 < slope2;

	// compare by y
	return line1_y < line2_y;
}


template<class t_hook, class t_vec, class t_line = std::pair<t_vec, t_vec>>
requires m::is_vec<t_vec>
struct IntersTreeLeaf
{
	using t_real = typename t_vec::value_type;

	const t_real *curX{nullptr};
	const std::vector<t_line> *lines{nullptr};
	std::size_t line_idx{0};

	t_real eps = std::numeric_limits<t_real>::epsilon();
	t_hook _h{};

	friend std::ostream& operator<<(std::ostream& ostr,
		const IntersTreeLeaf<t_hook, t_vec, t_line>& e)
	{
		ostr << std::get<0>((*e.lines)[e.line_idx])
			<< ", " << std::get<1>((*e.lines)[e.line_idx]);
		return ostr;
	}

	friend bool operator<(const IntersTreeLeaf<t_hook, t_vec, t_line>& e1,
		const IntersTreeLeaf<t_hook, t_vec, t_line>& e2)
	{
		const t_line& line1 = (*e1.lines)[e1.line_idx];
		const t_line& line2 = (*e2.lines)[e2.line_idx];

		return cmp_line<t_vec, t_line>(line1, line2, *e1.curX, e1.eps);
	}
};


/**
 * SO(2) rotation matrix
 */
template<class t_mat, class t_vec>
t_mat rotation_2d(const typename t_vec::value_type angle)
requires m::is_vec<t_vec> && m::is_mat<t_mat>
{
	using t_real = typename t_vec::value_type;

	const t_real c = std::cos(angle);
	const t_real s = std::sin(angle);

	return m::create<t_mat>({{c,s}, {-s,c}});
}


/**
 * line segment intersection via sweep
 * @returns [idx1, idx2, point]
 * @see (FUH 2020), ch. 2.3.2, pp. 69-80
 */
template<class t_vec, class t_line = std::tuple<t_vec, t_vec>,
	class t_real = typename t_vec::value_type>
requires m::is_vec<t_vec>
std::vector<std::tuple<std::size_t, std::size_t, t_vec>> intersect_sweep(
	const std::vector<t_line>& _lines, t_real eps = 1e-6)
{
	// additional parameter in t_line tuple serves as line group identifier
	constexpr bool use_line_groups = (std::tuple_size_v<t_line> > 2);
	using t_mat = m::mat<t_real, std::vector>;

	// look for vertical lines
	bool has_vert_line = 0;
	t_real min_angle_to_y{std::numeric_limits<t_real>::max()};
	std::vector<t_line> lines = _lines;

	for(const t_line& line : lines)
	{
		if(m::equals<t_real>(std::get<0>(line)[0], std::get<1>(line)[0], eps))
		{
			has_vert_line = 1;
		}
		else
		{
			// get angles relative to y axis
			t_real angle_to_y = line_angle<t_vec>(std::get<0>(line), std::get<1>(line)) + m::pi<t_real>/t_real(2);
			angle_to_y = m::mod_pos<t_real>(angle_to_y, t_real(2)*m::pi<t_real>);
			//std::cout << "angle: " << angle_to_y/m::pi<t_real>*180. << std::endl;

			if(angle_to_y > m::pi<t_real>/t_real(2))
				angle_to_y -= m::pi<t_real>;
			if(std::abs(angle_to_y) < std::abs(min_angle_to_y))
				min_angle_to_y = angle_to_y;
		}
	}

	// test rotation
	//has_vert_line = 1;
	//min_angle_to_y = m::pi<t_real> * 0.5;

	//if(has_vert_line)
	//	std::cout << "vertical line; next lowest angle: " << min_angle_to_y/m::pi<t_real>*180. << std::endl;

	// rotate all lines
	std::optional<t_mat> rotmat;
	if(has_vert_line)
	{
		rotmat = rotation_2d<t_mat, t_vec>(-min_angle_to_y * t_real(0.5));

		for(t_line& line : lines)
		{
			//std::cout << "line trafo: ";
			//print_line<t_vec, t_line>(std::cout, line);

			std::get<0>(line) = *rotmat * std::get<0>(line);
			std::get<1>(line) = *rotmat * std::get<1>(line);

			//std::cout << " -> ";
			//print_line<t_vec, t_line>(std::cout, line);
			//std::cout << std::endl;
		}
	}


	// order line vertices by x
	for(t_line& line : lines)
	{
		if(std::get<0>(line)[0] > std::get<1>(line)[0])
			std::swap(std::get<0>(line), std::get<1>(line));
	}


	enum class SweepEventType { LEFT_VERTEX, RIGHT_VERTEX, INTERSECTION };

	struct SweepEvent
	{
		t_real x;
		SweepEventType ty{SweepEventType::LEFT_VERTEX};

		std::size_t line_idx{};
		std::optional<std::size_t> lower_idx{}, upper_idx{};
		std::optional<t_vec> intersection{};

		void print(std::ostream& ostr)
		{
			std::string strty;
			if(ty == SweepEventType::LEFT_VERTEX)
				strty = "left_vertex";
			else if(ty == SweepEventType::RIGHT_VERTEX)
				strty = "right_vertex";
			else if(ty == SweepEventType::INTERSECTION)
				strty = "intersection";

			ostr << "x=" << std::setw(6) << x << ", type=" << std::setw(12)
				<< strty << ", line " << line_idx;
			if(lower_idx)
				ostr << ", lower=" << *lower_idx;
			if(upper_idx)
				ostr << ", upper=" << *upper_idx;

			if(intersection)
			{
				//using namespace m_ops;
				//ostr << ", intersection=" << *intersection;
				ostr << ", intersection=(" << (*intersection)[0] << "," << (*intersection)[1] << ")";
			}
		}
	};


	namespace intr = boost::intrusive;

	using t_leaf = IntersTreeLeaf<
		intr::avl_set_member_hook<
			intr::link_mode<
				intr::normal_link>>,
		t_vec, t_line>;

	using t_tree = intr::avltree<
		t_leaf, intr::member_hook<
			t_leaf, decltype(t_leaf::_h), &t_leaf::_h>>;


	// events
	auto events_comp = [](const SweepEvent& evt1, const SweepEvent& evt2) -> bool
	{
		return evt1.x > evt2.x;
	};

	std::priority_queue<SweepEvent, std::vector<SweepEvent>, decltype(events_comp)> events(events_comp);

	for(std::size_t line_idx=0; line_idx<lines.size(); ++line_idx)
	{
		const t_line& line = lines[line_idx];

		SweepEvent evtLeft{.x = std::get<0>(line)[0], .ty=SweepEventType::LEFT_VERTEX, .line_idx=line_idx};
		SweepEvent evtRight{.x = std::get<1>(line)[0], .ty=SweepEventType::RIGHT_VERTEX, .line_idx=line_idx};

		// wrong order of vertices?
		if(evtLeft.x > evtRight.x)
		{
			std::swap(evtLeft.ty, evtRight.ty);

			events.emplace(std::move(evtRight));
			events.emplace(std::move(evtLeft));
		}
		else
		{
			events.emplace(std::move(evtLeft));
			events.emplace(std::move(evtRight));
		}
	}


	auto add_intersection = [&events, &lines, eps]
		(std::size_t lower_idx, std::size_t upper_idx, t_real curX) -> void
		{
			//std::cout << "add_intersection: lower index: " << lower_idx << ", upper index: " << upper_idx << std::endl;
			const t_line& line1 = lines[lower_idx];
			const t_line& line2 = lines[upper_idx];

			if(auto [intersects, pt] = intersect_lines<t_line>(line1, line2, eps);
				intersects && !m::equals<t_real>(curX, pt[0], eps))
			{
				//std::cout << "add_intersection: intersection between lines " << lower_idx << " and " << upper_idx << std::endl;

				SweepEvent evtNext{.x = pt[0], .ty=SweepEventType::INTERSECTION,
					.lower_idx=lower_idx, .upper_idx=upper_idx,
					.intersection=pt};
				events.emplace(std::move(evtNext));
			}
		};

	// status
	t_tree status;

	// results
	std::vector<std::tuple<std::size_t, std::size_t, t_vec>> intersections;

	t_real curX = 0.;
	while(events.size())
	{
		SweepEvent evt{std::move(events.top())};
		events.pop();

		curX = evt.x;

		/*std::cout << "********* Event: "; evt.print(std::cout);
		std::cout << ", line order: ";
		for(auto theiter=status.begin(); theiter!=status.end(); theiter = std::next(theiter,1))
			std::cout << theiter->line_idx << ", ";
		std::cout << std::endl;*/

		switch(evt.ty)
		{
			/*
			 * arrived at left vertex
			 */
			case SweepEventType::LEFT_VERTEX:
			{
				// activate line
				t_leaf *leaf = new t_leaf{.curX=&curX, .lines=&lines, .line_idx=evt.line_idx, .eps=eps};
				auto iter = status.insert_equal(*leaf);

				auto iterPrev = (iter == status.begin() ? status.end() : std::prev(iter, 1));
				auto iterNext = (iter == status.end() ? status.end() : std::next(iter, 1));

				/*std::cout << "left vertex from line " << evt.line_idx;
				if(iterPrev != status.end()) std::cout << ", previous line: " << iterPrev->line_idx;
				if(iterNext != status.end()) std::cout << ", next line: " << iterNext->line_idx;
				std::cout << std::endl;*/

				// add possible intersection events
				if(iterPrev != iter && iterPrev != status.end())
					add_intersection(iterPrev->line_idx, evt.line_idx, curX);
				if(iterNext != iter && iterNext != status.end())
					add_intersection(evt.line_idx, iterNext->line_idx, curX);

				break;
			}

			/*
			 * arrived at right vertex
			 */
			case SweepEventType::RIGHT_VERTEX:
			{
				// find current line
				auto iter = std::find_if(status.begin(), status.end(),
					[&evt](const auto& leaf) -> bool
					{ return leaf.line_idx == evt.line_idx; });

				if(iter == status.end())
					continue;

				auto iterPrev = (iter == status.begin() ? status.end() : std::prev(iter, 1));
				auto iterNext = (iter == status.end() ? status.end() : std::next(iter, 1));

				// inactivate current line
				t_leaf *curLeaf = &*iter;
				iter = status.erase(iter);
				delete curLeaf;

				// add possible intersection event
				if(iterPrev != iterNext && iterPrev != status.end() && iterNext != status.end())
					add_intersection(iterPrev->line_idx, iterNext->line_idx, curX);

				break;
			}

			/*
			 * arrived at intersection
			 */
			case SweepEventType::INTERSECTION:
			{
				if(std::find_if(intersections.begin(), intersections.end(),
					[&evt, eps](const auto& inters) -> bool
						{ return m::equals<t_vec>(std::get<2>(inters), *evt.intersection, eps); })
					== intersections.end())
				{
					if constexpr(use_line_groups)
					{
						// if the lines belong to the same group, don't report the intersection
						if(std::get<2>(lines[*evt.lower_idx]) != std::get<2>(lines[*evt.upper_idx]))
						{
							// report an intersection
							intersections.emplace_back(std::make_tuple(*evt.lower_idx, *evt.upper_idx, *evt.intersection));
						}
					}
					else
					{
						// report an intersection
						intersections.emplace_back(std::make_tuple(*evt.lower_idx, *evt.upper_idx, *evt.intersection));
					}
				}
				else
				{
					// intersection already reported
					continue;
				}

				// find upper line
				auto iterUpper = std::find_if(status.begin(), status.end(),
					[&evt](const auto& leaf) -> bool
					{ return leaf.line_idx == evt.upper_idx; });

				// find lower line
				auto iterLower = std::find_if(status.begin(), status.end(),
					[&evt](const auto& leaf) -> bool
					{ return leaf.line_idx == evt.lower_idx; });

				if(iterUpper == status.end() || iterLower == status.end())
					continue;

				if(!cmp_line<t_vec, t_line>(lines[iterLower->line_idx], lines[iterUpper->line_idx], curX, eps))
				{
					std::swap(iterUpper->line_idx, iterLower->line_idx);
					std::swap(iterUpper, iterLower);
				}

				auto iterPrev = (iterUpper == status.begin() ? status.end() : std::prev(iterUpper, 1));
				auto iterNext = (iterLower == status.end() ? status.end() : std::next(iterLower, 1));

				// add possible intersection events
				if(iterPrev != iterUpper && iterPrev != status.end() && iterUpper != status.end())
					add_intersection(iterPrev->line_idx, iterUpper->line_idx, curX);
				if(iterNext != iterLower && iterNext != status.end() && iterLower != status.end())
					add_intersection(iterLower->line_idx, iterNext->line_idx, curX);

				break;
			}
		}


		// check leaf order
		/*auto lastiter = status.begin();
		for(auto theiter=status.begin(); theiter!=status.end(); theiter = std::next(theiter,1))
		{
			if(lastiter != theiter)
			{
				if(!cmp_line<t_vec, t_line>(lines[lastiter->line_idx],
					lines[theiter->line_idx], curX, eps))
				{
					std::cout << "Leaf order corrupted in lines " << lastiter->line_idx
						<< " and " << theiter->line_idx << ": ";

					const t_line& line1 = lines[lastiter->line_idx];
					const t_line& line2 = lines[theiter->line_idx];
					auto [slope1, offs1] = get_line_slope_offs<t_vec>(line1);
					auto [slope2, offs2] = get_line_slope_offs<t_vec>(line2);
					t_real line1_y = slope1*curX + offs1;
					t_real line2_y = slope2*curX + offs2;

					std::cout << "y1=" << line1_y << ", y2=" << line2_y
						<< ", slope1=" << slope1 << ", slope2=" << slope2 << std::endl;
				}
			}
			lastiter = theiter;
		}*/
	}


	// rotate intersection points back
	if(rotmat)
	{
		*rotmat = m::trans<t_mat>(*rotmat);

		for(auto& inters : intersections)
			std::get<2>(inters) = *rotmat * std::get<2>(inters);
	}

	//std::cout << std::endl;
	return intersections;
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// subvector
// ----------------------------------------------------------------------------

/**
 * maximum subvector
 * @see https://en.wikipedia.org/wiki/Maximum_subarray_problem
 */
template<class t_largernum, class t_arr>
std::tuple<std::size_t, std::size_t, t_largernum>
subvec_ineffic(const t_arr& arr)
requires m::is_iterable<t_arr> && m::is_basic_vec<t_arr>
{
	std::size_t start_idx = 0;
	std::size_t end_idx = 0;
	t_largernum val = -std::numeric_limits<t_largernum>::max();

	for(std::size_t start=0; start<arr.size(); ++start)
	{
		for(std::size_t end=start+1; end<=arr.size(); ++end)
		{
			t_largernum sum = std::accumulate(arr.begin()+start, arr.begin()+end, t_largernum{0});
			if(sum > val)
			{
				val = sum;
				start_idx = start;
				end_idx = end;
			}
		}
	}

	return std::make_tuple(start_idx, end_idx, val);
}


/**
 * maximum subvector
 * @see (FUH 2020), ch. 2.2.3, pp. 59-61
 * @see https://en.wikipedia.org/wiki/Maximum_subarray_problem
 */
template<class t_largernum, class t_arr>
std::tuple<std::size_t, std::size_t, t_largernum>
subvec_sweep(const t_arr& arr)
requires m::is_iterable<t_arr> && m::is_basic_vec<t_arr>
{
	std::size_t cached_start_idx = 0;
	std::size_t start_idx = 0;
	std::size_t end_idx = 0;

	t_largernum newval = 0;
	t_largernum val = 0;

	for(std::size_t idx=0; idx<arr.size(); ++idx)
	{
		if(newval < -arr[idx])
		{
			newval = t_largernum{0};
			cached_start_idx = idx+1;
		}
		else
		{
			newval += arr[idx];
		}

		if(newval > val)
		{
			val = newval;
			start_idx = cached_start_idx;
			end_idx = idx+1;
		}
	}

	return std::make_tuple(start_idx, end_idx, val);
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// closest pair
// @see (Klein 2005), ch. 2.2.2, pp. 53f; ch. 2.3.1, pp. 57f; ch. 2.4.1, pp. 93f
// @see (FUH 2020), ch. 2.2.2, pp. 58-69; ch. 2.3.1, pp. 62-69; ch. 2.4.1, pp. 95-96
// ----------------------------------------------------------------------------

template<class t_vec, class t_real = typename t_vec::value_type>
std::tuple<const t_vec*, const t_vec*, t_real>
closest_pair_ineff(const std::vector<t_vec>& points)
requires m::is_vec<t_vec>
{
	const t_vec* pt1 = nullptr;
	const t_vec* pt2 = nullptr;
	t_real dist = std::numeric_limits<t_real>::max();

	for(std::size_t i=0; i<points.size(); ++i)
	{
		for(std::size_t j=i+1; j<points.size(); ++j)
		{
			t_real norm = m::norm<t_vec>(points[i] - points[j]);
			if(norm < dist)
			{
				dist = norm;
				pt1 = &points[i];
				pt2 = &points[j];
			}
		}
	}

	return std::make_tuple(pt1, pt2, dist);
}


template<class t_hook, class T>
struct ClosestPairTreeLeaf
{
	const T* vec = nullptr;
	t_hook _h{};

	friend std::ostream& operator<<(std::ostream& ostr, const ClosestPairTreeLeaf<t_hook, T>& e)
	{
		ostr << *e.vec;
		return ostr;
	}

	friend bool operator<(const ClosestPairTreeLeaf<t_hook, T>& e1, const ClosestPairTreeLeaf<t_hook, T>& e2)
	{
		// compare by y
		return (*e1.vec)[1] < (*e2.vec)[1];
	}
};


/**
 * closest pair (2d)
 *  - @see (Klein 2005), ch 2.3.1, p. 57
 *  - @see (FUH 2020), ch. 2.3.1, pp. 62-69
 *	- @see https://en.wikipedia.org/wiki/Closest_pair_of_points_problem
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::tuple<t_vec, t_vec, t_real>
closest_pair_sweep(const std::vector<t_vec>& points)
requires m::is_vec<t_vec>
{
	namespace intr = boost::intrusive;

	using t_leaf = ClosestPairTreeLeaf<
		intr::avl_set_member_hook<
			intr::link_mode<
				intr::normal_link>>,
		t_vec>;
	using t_tree = intr::avltree<
		t_leaf, intr::member_hook<
			t_leaf, decltype(t_leaf::_h), &t_leaf::_h>>;


	std::vector<t_leaf> leaves;
	for(const t_vec& pt : points)
		leaves.emplace_back(t_leaf{.vec = &pt});

	std::stable_sort(leaves.begin(), leaves.end(),
	[](const t_leaf& leaf1, const t_leaf& leaf2) -> bool
	{
		// sort by x
		return (*leaf1.vec)[0] <= (*leaf2.vec)[0];
	});


	t_leaf* leaf1 = &leaves[0];
	t_leaf* leaf2 = &leaves[1];
	t_real dist = m::norm<t_vec>(*leaf1->vec - *leaf2->vec);

	t_tree status;
	status.insert_equal(*leaf1);
	status.insert_equal(*leaf2);


	auto iter1 = leaves.begin();
	for(auto iter2 = std::next(leaves.begin(), 2); iter2 != leaves.end(); )
	{
		if((*iter1->vec)[0] <= (*iter2->vec)[0]-dist)
		{
			// remove dead elements
			//std::cout << "Removing " << *iter1->vec << std::endl;
			status.erase(*iter1);
			std::advance(iter1, 1);
		}
		else
		{
			// insert newly active elements
			t_vec vec1 = *iter2->vec; vec1[1] -= dist;
			t_vec vec2 = *iter2->vec; vec2[1] += dist;
			auto [iterrange1, iterrange2] =
				status.bounded_range(t_leaf{.vec=&vec1}, t_leaf{.vec=&vec2}, 1, 1);

			for(auto iter=iterrange1; iter!=iterrange2; std::advance(iter,1))
			{
				t_real newdist = m::norm<t_vec>(*iter->vec - *iter2->vec);
				if(newdist < dist)
				{
					dist = newdist;
					leaf1 = &*iter;
					leaf2 = &*iter2;
				}
			}

			//std::cout << "Inserting " << *iter2->vec << std::endl;
			status.insert_equal(*iter2);
			std::advance(iter2, 1);
		}
	}

	return std::make_tuple(*leaf1->vec, *leaf2->vec, dist);
}


template<class t_vertex, class t_vec, std::size_t ...indices>
constexpr void _geo_vertex_assign(t_vertex& vert, const t_vec& vec, const std::index_sequence<indices...>&)
{
	namespace geo = boost::geometry;
	(geo::set<indices>(vert, vec[indices]), ...);
}


/**
 * closest pair (r-tree)
 */
template<std::size_t dim, class t_vec, class t_real = typename t_vec::value_type>
std::tuple<t_vec, t_vec, t_real>
closest_pair_rtree(const std::vector<t_vec>& _points)
requires m::is_vec<t_vec>
{
	std::vector<t_vec> points = _points;
	std::stable_sort(points.begin(), points.end(),
	[](const t_vec& pt1, const t_vec& pt2) -> bool
	{
		// sort by x
		return pt1[0] <= pt2[0];
	});


	namespace geo = boost::geometry;
	namespace geoidx = geo::index;

	using t_vertex = geo::model::point<t_real, dim, geo::cs::cartesian>;
	using t_box = geo::model::box<t_vertex>;
	using t_rtree_leaf = std::tuple<t_vertex, std::size_t>;
	using t_rtree = geoidx::rtree<t_rtree_leaf, geoidx::dynamic_linear>;

	t_rtree rt(typename t_rtree::parameters_type(points.size()));
	for(std::size_t ptidx=0; ptidx<points.size(); ++ptidx)
	{
		t_vertex vert;
		_geo_vertex_assign<t_vertex, t_vec>(vert, points[ptidx], std::make_index_sequence<dim>());

		rt.insert(std::make_tuple(vert, ptidx));
	}


	std::size_t idx1 = 0;
	std::size_t idx2 = 1;
	t_vec query1 = m::create<t_vec>(dim);
	t_vec query2 = m::create<t_vec>(dim);

	t_real dist = m::norm<t_vec>(points[idx2] - points[idx1]);
	for(std::size_t ptidx=1; ptidx<points.size(); ++ptidx)
	{
		query1[0] = points[ptidx][0] - dist;
		query2[0] = points[ptidx][0];

		for(std::size_t i=1; i<dim; ++i)
		{
			query1[i] = points[ptidx][i] - dist;
			query2[i] = points[ptidx][i] + dist;
		}

		t_vertex vert1, vert2;
		_geo_vertex_assign<t_vertex, t_vec>(vert1, query1, std::make_index_sequence<dim>());
		_geo_vertex_assign<t_vertex, t_vec>(vert2, query2, std::make_index_sequence<dim>());

		t_box query_obj(vert1, vert2);

		std::vector<t_rtree_leaf> query_answer;
		rt.query(geoidx::within(query_obj), std::back_inserter(query_answer));

		for(const auto& answ : query_answer)
		{
			std::size_t answidx = std::get<1>(answ);
			t_real newdist = m::norm<t_vec>(points[answidx] - points[ptidx]);
			if(newdist < dist)
			{
				dist = newdist;
				idx1 = answidx;
				idx2 = ptidx;
			}
		}
	}

	return std::make_tuple(points[idx1], points[idx2], dist);
}


/**
 * closest pair (range tree)
 * @see (FUH 2020), ch. 2.4.1, pp. 95-96; ch. 4.2.5, pp. 193-194
 */
template<std::size_t dim, class t_vec, class t_real = typename t_vec::value_type>
std::tuple<t_vec, t_vec, t_real>
closest_pair_rangetree(const std::vector<t_vec>& _points)
requires m::is_vec<t_vec>
{
	RangeTree<t_vec> tree;
	tree.insert(_points);

	// get x-sorted points
	std::vector<std::shared_ptr<const t_vec>> points;
	RangeTree<t_vec>::t_node::get_vecs(tree.get_root(), points);

	const t_vec* pt1 = points[0].get();
	const t_vec* pt2 = points[1].get();
	t_vec query1 = m::create<t_vec>(dim);
	t_vec query2 = m::create<t_vec>(dim);

	t_real dist = m::norm<t_vec>(*pt2 - *pt1);
	for(std::size_t ptidx=1; ptidx<points.size(); ++ptidx)
	{
		const t_vec& curpt = *points[ptidx];

		query1[0] = curpt[0] - dist;
		query2[0] = curpt[0];

		for(std::size_t i=1; i<dim; ++i)
		{
			query1[i] = curpt[i] - dist;
			query2[i] = curpt[i] + dist;
		}

		auto query_answer = tree.query_range(query1, query2);

		for(const auto& answ : query_answer)
		{
			t_real newdist = m::norm<t_vec>(*answ - curpt);
			if(answ.get() != &curpt && newdist < dist)
			{
				dist = newdist;
				pt1 = answ.get();
				pt2 = &curpt;
			}
		}
	}

	return std::make_tuple(*pt1, *pt2, dist);
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// trapezoid map
// Reference: (Berg 2008), ch. 6.2, pp. 128-133.
// ----------------------------------------------------------------------------

/**
 * trapezoid defined by a top and bottom line as well as a left and right point
 */
template<class t_vec> requires m::is_vec<t_vec>
class Trapezoid
{
public:
	using t_real = typename t_vec::value_type;
	using t_line = std::pair<t_vec, t_vec>;


public:
	Trapezoid() = default;
	~Trapezoid() = default;

	const t_vec& GetLeftPoint() const { return m_pointLeft; }
	const t_vec& GetRightPoint() const { return m_pointRight; }
	t_vec& GetLeftPoint() { return m_pointLeft; }
	t_vec& GetRightPoint() { return m_pointRight; }
	void SetLeftPoint(const t_vec& pt) { m_pointLeft = pt; }
	void SetRightPoint(const t_vec& pt) { m_pointRight = pt; }

	const t_line& GetTopLine() const { return m_lineTop; }
	const t_line& GetBottomLine() const { return m_lineBottom; }
	t_line& GetTopLine() { return m_lineTop; }
	t_line& GetBottomLine() { return m_lineBottom; }

	void SetTopLine(const t_line& line)
	{
		m_lineTop = line;

		// x component of line vertex 1 has to be left of vertex 2
		if(std::get<0>(m_lineTop)[0] > std::get<1>(m_lineTop)[0])
			std::swap(std::get<0>(m_lineTop), std::get<1>(m_lineTop));
	}

	void SetBottomLine(const t_line& line)
	{
		m_lineBottom = line;

		// x component of line vertex 1 has to be left of vertex 2
		if(std::get<0>(m_lineBottom)[0] > std::get<1>(m_lineBottom)[0])
			std::swap(std::get<0>(m_lineBottom), std::get<1>(m_lineBottom));
	}


	/**
	 * is point inside trapezoid (excluding border)?
	 */
	bool Contains(const t_vec& pt) const
	{
		// is point x left of left point x or right of right point x
		if(pt[0] <= m_pointLeft[0] || pt[0] >= m_pointRight[0])
			return false;

		// is point left of top line?
		if(side_of_line<t_vec, t_real>(
			std::get<0>(m_lineTop), std::get<1>(m_lineTop), pt) >= t_real(0))
			return false;

		// is point right of bottom line?
		if(side_of_line<t_vec, t_real>(
			std::get<0>(m_lineBottom), std::get<1>(m_lineBottom), pt) <= t_real(0))
			return false;

		return true;
	}


	/**
	 * trapezoid with no area
	 */
	bool IsEmpty(t_real eps=std::numeric_limits<t_real>::epsilon()) const
	{
		if(m::equals<t_real>(m_pointLeft[0], m_pointRight[0], eps))
			return true;
		if(is_line_equal<t_vec>(m_lineTop, m_lineBottom, eps))
			return true;

		return false;
	}


	/**
	 * let the trapezoid be the bounding box of the given points
	 */
	void SetBoundingBox(const std::vector<t_vec>& pts, t_real eps=std::numeric_limits<t_real>::epsilon())
	{
		t_real xmin = std::numeric_limits<t_real>::max();
		t_real xmax = -xmin;
		t_real ymin = xmin;
		t_real ymax = -ymin;

		for(const t_vec& pt : pts)
		{
			xmin = std::min(xmin, pt[0]);
			xmax = std::max(xmax, pt[0]);
			ymin = std::min(ymin, pt[1]);
			ymax = std::max(ymax, pt[1]);
		}

		xmin -= eps; xmax += eps;
		ymin -= eps; ymax += eps;

		m_pointLeft = m::create<t_vec>({xmin, ymin});
		m_pointRight = m::create<t_vec>({xmax, ymax});

		m_lineBottom = std::make_pair(m::create<t_vec>({xmin, ymin}), m::create<t_vec>({xmax, ymin}));
		m_lineTop = std::make_pair(m::create<t_vec>({xmin, ymax}), m::create<t_vec>({xmax, ymax}));
	}


	/**
	 * let the trapezoid be the bounding box of the given line segments
	 */
	void SetBoundingBox(const std::vector<t_line>& lines, t_real eps=std::numeric_limits<t_real>::epsilon())
	{
		std::vector<t_vec> pts;
		pts.reserve(lines.size()*2);

		for(const auto& line : lines)
		{
			pts.push_back(std::get<0>(line));
			pts.push_back(std::get<1>(line));
		}

		SetBoundingBox(pts, eps);
	}


private:
	t_vec m_pointLeft{}, m_pointRight{};
	t_line m_lineTop{}, m_lineBottom{};
};


enum class TrapezoidNodeType
{
	POINT,
	LINE,
	TRAPEZOID,
};


template<class t_vec> requires m::is_vec<t_vec>
class TrapezoidNode
{
public:
	TrapezoidNode() = default;
	virtual ~TrapezoidNode() = default;

	virtual TrapezoidNodeType GetType() const = 0;
	virtual bool IsLeft(const t_vec& vec) const = 0;

	virtual const std::shared_ptr<TrapezoidNode<t_vec>> GetLeft() const = 0;
	virtual const std::shared_ptr<TrapezoidNode<t_vec>> GetRight() const = 0;

	virtual void SetLeft(const std::shared_ptr<TrapezoidNode<t_vec>>& left) = 0;
	virtual void SetRight(const std::shared_ptr<TrapezoidNode<t_vec>>& right) = 0;
};


template<class t_vec> requires m::is_vec<t_vec>
class TrapezoidNodePoint : public TrapezoidNode<t_vec>
{
public:
	TrapezoidNodePoint(const t_vec& vec = m::zero<t_vec>()) : m_vec{vec}
	{}

	TrapezoidNodePoint(
		const std::shared_ptr<TrapezoidNode<t_vec>>& left,
		const std::shared_ptr<TrapezoidNode<t_vec>>& right)
		: m_left{left}, m_right{right}
	{}

	virtual ~TrapezoidNodePoint() = default;

	virtual const std::shared_ptr<TrapezoidNode<t_vec>> GetLeft() const override
	{
		return m_left;
	}

	virtual const std::shared_ptr<TrapezoidNode<t_vec>> GetRight() const override
	{
		return m_right;
	}

	virtual void SetLeft(const std::shared_ptr<TrapezoidNode<t_vec>>& left) override
	{
		m_left = left;
	}

	virtual void SetRight(const std::shared_ptr<TrapezoidNode<t_vec>>& right) override
	{
		m_right = right;
	}

	virtual TrapezoidNodeType GetType() const override
	{
		return TrapezoidNodeType::POINT;
	}

	virtual bool IsLeft(const t_vec& vec) const override
	{
		// is vec left of m_vec?
		return vec[0] <= m_vec[0];
	}


	const t_vec& GetPoint() const { return m_vec; }
	t_vec& GetPoint() { return m_vec; }
	void SetPoint(const t_vec& vec) { m_vec = vec; }


private:
	std::shared_ptr<TrapezoidNode<t_vec>> m_left{}, m_right{};
	t_vec m_vec{};
};


template<class t_vec> requires m::is_vec<t_vec>
class TrapezoidNodeLine : public TrapezoidNode<t_vec>
{
public:
	using t_real = typename t_vec::value_type;
	using t_line = std::pair<t_vec, t_vec>;


public:
	TrapezoidNodeLine(const t_line& line
		= std::make_pair<t_vec, t_vec>(m::zero<t_vec>(), m::zero<t_vec>()))
		: m_line{line}
	{}

	TrapezoidNodeLine(
		const std::shared_ptr<TrapezoidNode<t_vec>>& left,
		const std::shared_ptr<TrapezoidNode<t_vec>>& right)
		: m_left{left}, m_right{right}
	{}

	virtual ~TrapezoidNodeLine() = default;


	virtual const std::shared_ptr<TrapezoidNode<t_vec>> GetLeft() const override
	{
		return m_left;
	}

	virtual const std::shared_ptr<TrapezoidNode<t_vec>> GetRight() const override
	{
		return m_right;
	}

	virtual void SetLeft(const std::shared_ptr<TrapezoidNode<t_vec>>& left) override
	{
		m_left = left;
	}

	virtual void SetRight(const std::shared_ptr<TrapezoidNode<t_vec>>& right) override
	{
		m_right = right;
	}

	virtual TrapezoidNodeType GetType() const override
	{
		return TrapezoidNodeType::LINE;
	}

	virtual bool IsLeft(const t_vec& vec) const override
	{
		// is vec left of m_line?
		return side_of_line<t_vec, t_real>(
			std::get<0>(m_line), std::get<1>(m_line), vec) >= t_real(0);
	}


	const t_line& GetLine() const { return m_line; }
	t_line& GetLine() { return m_line; }
	void SetLine(const t_line& line) { m_line = line; }


private:
	std::shared_ptr<TrapezoidNode<t_vec>> m_left{}, m_right{};
	t_line m_line{};
};


template<class t_vec> requires m::is_vec<t_vec>
class TrapezoidNodeTrapezoid : public TrapezoidNode<t_vec>
{
public:
	TrapezoidNodeTrapezoid(const std::shared_ptr<Trapezoid<t_vec>>& trapezoid = nullptr)
		: m_trapezoid{trapezoid}
	{
	}

	virtual ~TrapezoidNodeTrapezoid() = default;

	virtual TrapezoidNodeType GetType() const override
	{
		return TrapezoidNodeType::TRAPEZOID;
	}


	const std::shared_ptr<Trapezoid<t_vec>> GetTrapezoid() const
	{
		return m_trapezoid;
	}

	void SetTrapezoid(const std::shared_ptr<Trapezoid<t_vec>>& trapezoid)
	{
		return m_trapezoid = trapezoid;
	}


protected:
	virtual const std::shared_ptr<TrapezoidNode<t_vec>> GetLeft() const override
	{ return nullptr; }

	virtual const std::shared_ptr<TrapezoidNode<t_vec>> GetRight() const override
	{ return nullptr; }

	virtual void SetLeft(const std::shared_ptr<TrapezoidNode<t_vec>>&) override
	{ }

	virtual void SetRight(const std::shared_ptr<TrapezoidNode<t_vec>>&) override
	{ }

	virtual bool IsLeft(const t_vec&) const override
	{ return false; }


private:
	std::shared_ptr<Trapezoid<t_vec>> m_trapezoid{};
};


template<class t_vec, class t_line = std::pair<t_vec, t_vec>>
requires m::is_vec<t_vec>
std::ostream& operator<<(std::ostream& ostr,
	const std::pair<std::shared_ptr<TrapezoidNode<t_vec>>, int>& node_depth)
{
	const auto& node = *std::get<0>(node_depth);
	int depth = std::get<1>(node_depth);

	auto print_indent = [&ostr, depth]() -> void
	{
		for(int i=0; i<depth; ++i)
			ostr << "  ";
	};

	print_indent();

	switch(node.GetType())
	{
		case TrapezoidNodeType::POINT:
		{
			auto ptnode = dynamic_cast<const TrapezoidNodePoint<t_vec>&>(node);
			const auto& pt = ptnode.GetPoint();
			ostr << "point: ";
			print_point<t_vec>(ostr, pt);

			break;
		}

		case TrapezoidNodeType::LINE:
		{
			auto linenode = dynamic_cast<const TrapezoidNodeLine<t_vec>&>(node);
			const auto& line = linenode.GetLine();
			ostr << "line: ";
			print_line<t_vec>(ostr, line);

			break;
		}

		case TrapezoidNodeType::TRAPEZOID:
		{
			auto trnode = dynamic_cast<const TrapezoidNodeTrapezoid<t_vec>&>(node);

			ostr << "trapezoid: ";
			ostr << "left: ";
			print_point<t_vec>(ostr, trnode.GetTrapezoid()->GetLeftPoint());
			ostr << ", right: ";
			print_point<t_vec>(ostr, trnode.GetTrapezoid()->GetRightPoint());
			ostr << ", bottom: ";
			print_line<t_vec>(ostr, trnode.GetTrapezoid()->GetBottomLine());
			ostr << ", top: ";
			print_line<t_vec>(ostr, trnode.GetTrapezoid()->GetTopLine());

			break;
		}
	};

	ostr << "\n";

	if(node.GetLeft())
	{
		print_indent();
		ostr << "left node:\n";
		ostr << std::make_pair(node.GetLeft(), depth+1);
	}
	if(node.GetRight())
	{
		print_indent();
		ostr << "right node:\n";
		ostr << std::make_pair(node.GetRight(), depth+1);
	}

	return ostr;
}


/**
 * find a neighbouring trapezoid in the tree
 */
template<class t_vec, class t_line=std::pair<t_vec, t_vec>, class t_real=typename t_vec::value_type>
requires m::is_vec<t_vec>
std::shared_ptr<TrapezoidNodeTrapezoid<t_vec>> find_neighbour_trapezoid(
	const std::shared_ptr<TrapezoidNode<t_vec>>& node,
	const std::shared_ptr<Trapezoid<t_vec>>& trap,
	bool left = 1, bool top = 1,
	t_real eps=std::numeric_limits<t_real>::epsilon())
{
	const t_line* lineTop = &trap->GetTopLine();
	const t_line* lineBottom = &trap->GetBottomLine();
	const t_vec* ptLeft = &trap->GetLeftPoint();
	const t_vec* ptRight = &trap->GetRightPoint();

	// possible neighbours
	std::vector<std::shared_ptr<TrapezoidNodeTrapezoid<t_vec>>> candidates;

	// function to traverse the tree
	std::function<void(const std::shared_ptr<TrapezoidNode<t_vec>>&)> traverse;
	traverse = [&candidates, &traverse, lineTop, lineBottom, ptLeft, ptRight, left, top, eps]
		(const std::shared_ptr<TrapezoidNode<t_vec>>& node) -> void
	{
		// TODO: better use of binary search tree, no need to check all nodes...
		if(node->GetLeft())
			traverse(node->GetLeft());
		if(node->GetRight())
			traverse(node->GetRight());

		if(node->GetType() == TrapezoidNodeType::TRAPEZOID)
		{
			auto trnode = std::dynamic_pointer_cast<TrapezoidNodeTrapezoid<t_vec>>(node);
			const t_line& lineTop2 = trnode->GetTrapezoid()->GetTopLine();
			const t_line& lineBottom2 = trnode->GetTrapezoid()->GetBottomLine();
			const t_vec& ptLeft2 = trnode->GetTrapezoid()->GetLeftPoint();
			const t_vec& ptRight2 = trnode->GetTrapezoid()->GetRightPoint();

			if(left && top)
			{
				if(is_line_equal<t_vec>(*lineTop, lineTop2, eps)
					&& m::equals<t_vec>(*ptLeft, ptRight2, eps))
					candidates.push_back(trnode);
			}
			else if(left && !top)
			{
				if(is_line_equal<t_vec>(*lineBottom, lineBottom2, eps)
					&& m::equals<t_vec>(*ptLeft, ptRight2, eps))
					candidates.push_back(trnode);
			}
			else if(!left && top)
			{
				if(is_line_equal<t_vec>(*lineTop, lineTop2, eps)
					&& m::equals<t_vec>(*ptRight, ptLeft2, eps))
					candidates.push_back(trnode);
			}
			else if(!left && !top)
			{
				if(is_line_equal<t_vec>(*lineBottom, lineBottom2, eps)
					&& m::equals<t_vec>(*ptRight, ptLeft2, eps))
					candidates.push_back(trnode);
			}
		}
	};

	// search the tree
	traverse(node);

	// no neighbours found
	if(candidates.size() == 0)
		return nullptr;

	if(candidates.size() > 1)
	{
		std::cerr << __PRETTY_FUNCTION__
			<< ": Warning: more than one neighbour found." << std::endl;
	}

	return candidates[0];
}


/**
 * find trapezoid containing a given point
 */
template<class t_vec, class t_real=typename t_vec::value_type>
requires m::is_vec<t_vec>
std::shared_ptr<TrapezoidNodeTrapezoid<t_vec>> find_trapezoid(
	const std::shared_ptr<TrapezoidNode<t_vec>>& node, const t_vec& pt)
{
	// possible trapezoids
	std::vector<std::shared_ptr<TrapezoidNodeTrapezoid<t_vec>>> candidates;

	// function to traverse the tree
	std::function<void(const std::shared_ptr<TrapezoidNode<t_vec>>&)> traverse;
	traverse = [&candidates, &pt, &traverse]
		(const std::shared_ptr<TrapezoidNode<t_vec>>& node) -> void
	{
		if(node->IsLeft(pt))
		{
			if(node->GetLeft())
				traverse(node->GetLeft());
		}
		else
		{
			if(node->GetRight())
				traverse(node->GetRight());
		}

		if(node->GetType() == TrapezoidNodeType::TRAPEZOID)
		{
			auto trnode = std::dynamic_pointer_cast<TrapezoidNodeTrapezoid<t_vec>>(node);
			if(trnode->GetTrapezoid()->Contains(pt))
				candidates.push_back(trnode);
		}
	};

	// search the tree
	traverse(node);

	// no neighbours found
	if(candidates.size() == 0)
		return nullptr;

	if(candidates.size() > 1)
	{
		std::cerr << __PRETTY_FUNCTION__
			<< ": Warning: more than one trapezoid found." << std::endl;
	}

	return candidates[0];
}


/**
 * replace old trapezoid node pointers -> new node pointers
 */
template<class t_vec> requires m::is_vec<t_vec>
void replace_trapezoid_ptr(std::shared_ptr<TrapezoidNode<t_vec>> node,
	std::shared_ptr<TrapezoidNodeTrapezoid<t_vec>> old_node,
	std::shared_ptr<TrapezoidNode<t_vec>> new_node)
{
	if(!node || !old_node)
		return;

	auto old_trap = old_node->GetTrapezoid();
	if(!old_trap)
		return;

	if(node->GetLeft() && node->GetLeft()->GetType() == TrapezoidNodeType::TRAPEZOID)
	{
		auto trnode = std::dynamic_pointer_cast<TrapezoidNodeTrapezoid<t_vec>>(node->GetLeft());

		if(trnode->GetTrapezoid() == old_trap)
			node->SetLeft(new_node);
	}
	else
	{
		replace_trapezoid_ptr<t_vec>(node->GetLeft(), old_node, new_node);
	}

	if(node->GetRight() && node->GetRight()->GetType() == TrapezoidNodeType::TRAPEZOID)
	{
		auto trnode = std::dynamic_pointer_cast<TrapezoidNodeTrapezoid<t_vec>>(node->GetRight());

		if(trnode->GetTrapezoid() == old_trap)
			node->SetRight(new_node);
	}
	else
	{
		replace_trapezoid_ptr<t_vec>(node->GetRight(), old_node, new_node);
	}
}


/**
 * cut the bottom and top lines to not exceed the left and right point
 */
template<class t_vec, class t_line = std::pair<t_vec, t_vec>, class t_real = typename t_vec::value_type>
requires m::is_vec<t_vec>
void fit_trapezoid_lines(std::shared_ptr<TrapezoidNode<t_vec>> node)
{
	if(node->GetLeft())
		fit_trapezoid_lines<t_vec>(node->GetLeft());
	if(node->GetRight())
		fit_trapezoid_lines<t_vec>(node->GetRight());

	if(node->GetType() == TrapezoidNodeType::TRAPEZOID)
	{
		auto trnode = std::dynamic_pointer_cast<TrapezoidNodeTrapezoid<t_vec>>(node);
		auto trap = trnode->GetTrapezoid();

		t_real x0 = trap->GetLeftPoint()[0];
		t_real x1 = trap->GetRightPoint()[0];

		t_line lineBottom = trap->GetBottomLine();
		t_real bottom_y0 = get_line_y<t_vec>(lineBottom, x0);
		t_real bottom_y1 = get_line_y<t_vec>(lineBottom, x1);
		std::get<0>(lineBottom)[0] = x0;
		std::get<0>(lineBottom)[1] = bottom_y0;
		std::get<1>(lineBottom)[0] = x1;
		std::get<1>(lineBottom)[1] = bottom_y1;;
		trap->SetBottomLine(lineBottom);

		t_line lineTop = trap->GetTopLine();
		t_real top_y0 = get_line_y<t_vec>(lineTop, x0);
		t_real top_y1 = get_line_y<t_vec>(lineTop, x1);
		std::get<0>(lineTop)[0] = x0;
		std::get<0>(lineTop)[1] = top_y0;
		std::get<1>(lineTop)[0] = x1;
		std::get<1>(lineTop)[1] = top_y1;
		trap->SetTopLine(lineTop);
	}
}


/**
 * transform all points, lines and trapezoids in the tree
 */
template<class t_vec, class t_mat>
requires m::is_vec<t_vec> && m::is_mat<t_mat>
void trafo_trapezoid_tree(std::shared_ptr<TrapezoidNode<t_vec>> node, const t_mat& mat,
	std::shared_ptr<std::unordered_set<void*>> cache=nullptr)
{
	// prevent pointers in the tree to be transformed multiple times
	if(!cache)
		cache = std::make_shared<std::unordered_set<void*>>();

	if(node->GetType() == TrapezoidNodeType::POINT)
	{
		auto ptnode = std::dynamic_pointer_cast<TrapezoidNodePoint<t_vec>>(node);

		if(cache->find(ptnode.get()) == cache->end())
		{
			ptnode->SetPoint(mat * ptnode->GetPoint());
			cache->insert(ptnode.get());
		}
	}
	else if(node->GetType() == TrapezoidNodeType::LINE)
	{
		auto lnnode = std::dynamic_pointer_cast<TrapezoidNodeLine<t_vec>>(node);

		if(cache->find(lnnode.get()) == cache->end())
		{
			std::get<0>(lnnode->GetLine()) = mat * std::get<0>(lnnode->GetLine());
			std::get<1>(lnnode->GetLine()) = mat * std::get<1>(lnnode->GetLine());

			cache->insert(lnnode.get());
		}
	}
	else if(node->GetType() == TrapezoidNodeType::TRAPEZOID)
	{
		auto trnode = std::dynamic_pointer_cast<TrapezoidNodeTrapezoid<t_vec>>(node);
		auto trap = trnode->GetTrapezoid();

		if(cache->find(trap.get()) == cache->end())
		{
			trap->SetLeftPoint(mat * trap->GetLeftPoint());
			trap->SetRightPoint(mat * trap->GetRightPoint());
			std::get<0>(trap->GetTopLine()) = mat * std::get<0>(trap->GetTopLine());
			std::get<1>(trap->GetTopLine()) = mat * std::get<1>(trap->GetTopLine());
			std::get<0>(trap->GetBottomLine()) = mat * std::get<0>(trap->GetBottomLine());
			std::get<1>(trap->GetBottomLine()) = mat * std::get<1>(trap->GetBottomLine());

			cache->insert(trap.get());
		}
	}

	if(node->GetLeft())
		trafo_trapezoid_tree<t_vec, t_mat>(node->GetLeft(), mat, cache);
	if(node->GetRight())
		trafo_trapezoid_tree<t_vec, t_mat>(node->GetRight(), mat, cache);
}


/**
 * save the trapezoid tree as an svg
 */
template<class t_vec, class t_line = std::pair<t_vec, t_vec>, class t_real = typename t_vec::value_type>
requires m::is_vec<t_vec>
void save_trapezoid_svg(const std::shared_ptr<TrapezoidNode<t_vec>>& node,
	const std::string& file, const std::vector<t_line>* lines = nullptr)
{
	namespace geo = boost::geometry;
	using t_geovertex = geo::model::point<t_real, 2, geo::cs::cartesian>;
	using t_geopoly = geo::model::polygon<t_geovertex, true /*cw*/, false /*closed*/>;
	using t_geoline = geo::model::linestring<t_geovertex>;
	using t_geosvg = geo::svg_mapper<t_geovertex>;

	// the tree is a dag -> avoid writing the same pointer several times
	std::unordered_set<void*> cache;

	std::ofstream ofstr{file};
	t_geosvg svg{ofstr, 100, 100, "width = \"500px\" height = \"500px\""};

	// function to traverse the tree
	std::function<void(const std::shared_ptr<TrapezoidNode<t_vec>>&)> traverse;
	traverse = [&svg, &traverse, &cache]
	(const std::shared_ptr<TrapezoidNode<t_vec>>& node) -> void
	{
		if(node->GetLeft())
			traverse(node->GetLeft());
		if(node->GetRight())
			traverse(node->GetRight());

		if(node->GetType() == TrapezoidNodeType::TRAPEZOID)
		{
			auto trnode = std::dynamic_pointer_cast<TrapezoidNodeTrapezoid<t_vec>>(node);
			auto trap = trnode->GetTrapezoid();
			if(cache.find(trap.get()) == cache.end())
			{
				cache.insert(trap.get());

				auto lineTop = trap->GetTopLine();
				auto lineBottom = trap->GetBottomLine();

				t_geopoly poly;
				poly.outer().push_back(t_geovertex{std::get<0>(lineTop)[0], std::get<0>(lineTop)[1]});
				poly.outer().push_back(t_geovertex{std::get<0>(lineBottom)[0], std::get<0>(lineBottom)[1]});
				poly.outer().push_back(t_geovertex{std::get<1>(lineBottom)[0], std::get<1>(lineBottom)[1]});
				poly.outer().push_back(t_geovertex{std::get<1>(lineTop)[0], std::get<1>(lineTop)[1]});

				svg.add(poly);
				svg.map(poly, "stroke: #000000; stroke-width: 1px; fill: none;", 1.);
			}
		}
	};

	traverse(node);

	if(lines)
	{
		for(const auto& line : *lines)
		{
			t_geoline theline;
			theline.push_back(t_geovertex{std::get<0>(line)[0], std::get<0>(line)[1]});
			theline.push_back(t_geovertex{std::get<1>(line)[0], std::get<1>(line)[1]});

			svg.add(theline);
			svg.map(theline, "stroke: #ff0000; stroke-width: 2px; fill: none;", 1.);
		}
	}

	ofstr.flush();
}


/**
 * check if the line segments share the same x coordinates
 */
template<class t_vec, class t_line = std::pair<t_vec, t_vec>,
	class t_real = typename t_vec::value_type>
requires m::is_vec<t_vec>
bool check_line_equal_x(const std::vector<t_line>& lines,
	bool exclude_duplicates = true,  // exclude coinciding points
	t_real eps=std::numeric_limits<t_real>::epsilon())
{
	struct CompareReals
	{
		CompareReals() = delete;
		CompareReals(t_real eps) : m_eps{eps}
		{}

		bool operator()(t_real x, t_real y) const
		{
			//std::cout << x << " < " << y << ", eps = " << m_eps << std::endl;

			if(m::equals<t_real>(x, y, m_eps))
				return false;

			return x < y;
		}

		t_real m_eps = std::numeric_limits<t_real>::epsilon();
	};

	std::vector<t_vec> alreadyseen0, alreadyseen1;
	std::set<t_real, CompareReals> cache_x0{CompareReals{eps}},
		cache_x1{CompareReals{eps}};

	for(const t_line& line : lines)
	{
		const t_vec& pt0 = std::get<0>(line);
		const t_vec& pt1 = std::get<1>(line);

		if(exclude_duplicates)
		{
			if(std::find_if(alreadyseen0.begin(), alreadyseen0.end(),
				[&pt0, eps](const t_vec& pt) -> bool
				{
					return m::equals<t_vec>(pt0, pt, eps);
				}) != alreadyseen0.end())
				continue;

			if(std::find_if(alreadyseen1.begin(), alreadyseen1.end(),
				[&pt1, eps](const t_vec& pt) -> bool
				{
					return m::equals<t_vec>(pt1, pt, eps);
				}) != alreadyseen1.end())
				continue;
		}

		t_real x0 = pt0[0];
		t_real x1 = pt1[0];

		if(cache_x0.find(x0) != cache_x0.end())
			return true;
		if(cache_x1.find(x1) != cache_x1.end())
			return true;

		cache_x0.insert(x0);
		cache_x1.insert(x1);

		if(exclude_duplicates)
		{
			alreadyseen0.push_back(pt0);
			alreadyseen1.push_back(pt1);
		}
	}

	return false;
}


/**
 * try to unite two adjacent trapezoids
 */
template<class t_vec, class t_real = typename t_vec::value_type>
requires m::is_vec<t_vec>
bool try_unite_trapezoids(
	std::shared_ptr<Trapezoid<t_vec>>* trap1, // left trapezoid
	std::shared_ptr<Trapezoid<t_vec>>* trap2, // right trapezoid
	t_real eps=std::numeric_limits<t_real>::epsilon())
{
	if(!*trap1 || !*trap2)
		return false;
	if(*trap1 == *trap2)
		return false;

	if((*trap1)->GetLeftPoint()[0] > (*trap2)->GetLeftPoint()[0])
		std::swap(trap1, trap2);

	if(!is_line_equal<t_vec>((*trap1)->GetTopLine(), (*trap2)->GetTopLine(), eps))
		return false;
	if(!is_line_equal<t_vec>((*trap1)->GetBottomLine(), (*trap2)->GetBottomLine(), eps))
		return false;

	if(!m::equals<t_real>((*trap1)->GetRightPoint()[0], (*trap2)->GetLeftPoint()[0], eps))
		return false;
	//if(!m::equals<t_real>((*trap1)->GetRightPoint()[1], (*trap2)->GetLeftPoint()[1], eps))
	//	return false;

	(*trap1)->SetRightPoint((*trap2)->GetRightPoint());
	*trap2 = *trap1;

	return true;
}


/**
 * remove empty nodes
 * @return can node be deleted?
 */
template<class t_vec, class t_real = typename t_vec::value_type>
requires m::is_vec<t_vec>
bool clean_trapezoid_tree(std::shared_ptr<TrapezoidNode<t_vec>> node)
{
	if(!node)
		return true;

	// don't remove trapezoid nodes
	if(node->GetType() == TrapezoidNodeType::TRAPEZOID)
		return false;

	// delete point or line nodes with no children
	if(!node->GetLeft() && !node->GetRight())
		return true;

	if(clean_trapezoid_tree<t_vec>(node->GetLeft()))
		node->SetLeft(nullptr);
	if(clean_trapezoid_tree<t_vec>(node->GetRight()))
		node->SetRight(nullptr);

	return !node->GetLeft() && !node->GetRight();
}


/**
 * create a trapezoid tree
 * @see (Berg 2008), pp. 128-133 and pp. 137-139
 */
template<class t_vec, class t_line = std::pair<t_vec, t_vec>, class t_real = typename t_vec::value_type>
requires m::is_vec<t_vec>
std::shared_ptr<TrapezoidNode<t_vec>>
create_trapezoid_tree(const std::vector<t_line>& _lines,
	bool randomise=true, bool shear=true,
	t_real padding=1., t_real eps=1e-5)
{
	using t_mat = m::mat<t_real, std::vector>;
	std::vector<t_line> lines = _lines;

	t_real shear_eps = eps * 10.;
	if(shear)
	{
		// shear until no more x coordinates coincide
		t_real shear_mult = 1;
		while(check_line_equal_x<t_vec>(lines, true, eps))
		{
			t_mat shear = m::shear<t_mat>(2, 2, 0, 1, shear_eps);

			for(auto& line : lines)
			{
				std::get<0>(line) = shear*std::get<0>(line);
				std::get<1>(line) = shear*std::get<1>(line);
			}

			++shear_mult;
		}
		shear_eps *= shear_mult;
	}

	auto box = std::make_shared<Trapezoid<t_vec>>();
	box->SetBoundingBox(lines, padding);

	std::shared_ptr<TrapezoidNode<t_vec>> root =
		std::make_shared<TrapezoidNodeTrapezoid<t_vec>>(box);

	if(randomise)
		std::ranges::shuffle(lines, std::mt19937{std::random_device{}()});

	// add lines to the tree
	for(std::size_t lineidx=0; lineidx<lines.size(); ++lineidx)
	{
		const auto& line = lines[lineidx];
		std::vector<std::shared_ptr<TrapezoidNodeTrapezoid<t_vec>>> intersecting_trapezoids;

		const t_vec& leftpt = std::get<0>(line);
		const t_vec& rightpt = std::get<1>(line);

		if(auto trap_node = find_trapezoid<t_vec>(root, leftpt); trap_node)
			intersecting_trapezoids.push_back(trap_node);

		if(intersecting_trapezoids.size() == 0)
			continue;

		auto cur_trap = intersecting_trapezoids[0];
		while(cur_trap)
		{
			const auto& trap_rightpt = cur_trap->GetTrapezoid()->GetRightPoint();
			if(rightpt[0] <= trap_rightpt[0])
				break;

			if(side_of_line<t_vec, t_real>(leftpt, rightpt, trap_rightpt) >= t_real(0))
				cur_trap = find_neighbour_trapezoid<t_vec>(root, cur_trap->GetTrapezoid(), 0, 0, eps);
			else
				cur_trap = find_neighbour_trapezoid<t_vec>(root, cur_trap->GetTrapezoid(), 0, 1, eps);

			if(cur_trap)
				intersecting_trapezoids.push_back(cur_trap);
		}


		auto create_trapnode = [](std::shared_ptr<Trapezoid<t_vec>> trap)
			-> std::shared_ptr<TrapezoidNodeTrapezoid<t_vec>>
		{
			return std::make_shared<TrapezoidNodeTrapezoid<t_vec>>(trap);
		};

		if(intersecting_trapezoids.size() == 0)
		{
			break;
		}
		else if(intersecting_trapezoids.size() == 1)
		{
			auto cur_trap_node = intersecting_trapezoids[0];
			auto cur_trap = cur_trap_node->GetTrapezoid();

			auto trap_left = std::make_shared<Trapezoid<t_vec>>();
			trap_left->SetLeftPoint(cur_trap->GetLeftPoint());
			trap_left->SetRightPoint(leftpt);
			trap_left->SetTopLine(cur_trap->GetTopLine());
			trap_left->SetBottomLine(cur_trap->GetBottomLine());

			auto trap_right = std::make_shared<Trapezoid<t_vec>>();
			trap_right->SetLeftPoint(rightpt);
			trap_right->SetRightPoint(cur_trap->GetRightPoint());
			trap_right->SetTopLine(cur_trap->GetTopLine());
			trap_right->SetBottomLine(cur_trap->GetBottomLine());

			auto trap_top = std::make_shared<Trapezoid<t_vec>>();
			trap_top->SetLeftPoint(leftpt);
			trap_top->SetRightPoint(rightpt);
			trap_top->SetTopLine(cur_trap->GetTopLine());
			trap_top->SetBottomLine(line);

			auto trap_bottom = std::make_shared<Trapezoid<t_vec>>();
			trap_bottom->SetLeftPoint(leftpt);
			trap_bottom->SetRightPoint(rightpt);
			trap_bottom->SetTopLine(line);
			trap_bottom->SetBottomLine(cur_trap->GetBottomLine());

			auto trap_left_node = create_trapnode(trap_left);
			auto trap_right_node = create_trapnode(trap_right);
			auto trap_top_node = create_trapnode(trap_top);
			auto trap_bottom_node = create_trapnode(trap_bottom);

			auto line_node = std::make_shared<TrapezoidNodeLine<t_vec>>(line);
			if(!trap_top->IsEmpty(eps))
				line_node->SetLeft(trap_top_node);
			if(!trap_bottom->IsEmpty(eps))
				line_node->SetRight(trap_bottom_node);

			auto rightpt_node = std::make_shared<TrapezoidNodePoint<t_vec>>(rightpt);
			rightpt_node->SetLeft(line_node);
			if(!trap_right->IsEmpty(eps))
			rightpt_node->SetRight(trap_right_node);

			auto leftpt_node = std::make_shared<TrapezoidNodePoint<t_vec>>(leftpt);
			leftpt_node->SetLeft(trap_left_node);
			leftpt_node->SetRight(rightpt_node);

			fit_trapezoid_lines<t_vec>(leftpt_node);

			if(root->GetType() == TrapezoidNodeType::TRAPEZOID &&
				cur_trap == std::dynamic_pointer_cast<TrapezoidNodeTrapezoid<t_vec>>(root)->GetTrapezoid())
			{
				// replace root node
				root = leftpt_node;
			}
			else
			{
				// replace child node pointers
				replace_trapezoid_ptr<t_vec>(root, cur_trap_node, leftpt_node);
			}
		}
		else if(intersecting_trapezoids.size() > 1)
		{
			// first trapezoid
			auto first_trap_node = *intersecting_trapezoids.begin();
			auto first_trap = first_trap_node->GetTrapezoid();

			auto first_left = std::make_shared<Trapezoid<t_vec>>();
			first_left->SetLeftPoint(first_trap->GetLeftPoint());
			first_left->SetRightPoint(leftpt);
			first_left->SetTopLine(first_trap->GetTopLine());
			first_left->SetBottomLine(first_trap->GetBottomLine());

			auto first_top = std::make_shared<Trapezoid<t_vec>>();
			first_top->SetLeftPoint(leftpt);
			first_top->SetRightPoint(first_trap->GetRightPoint());
			first_top->SetTopLine(first_trap->GetTopLine());
			first_top->SetBottomLine(line);

			auto first_bottom = std::make_shared<Trapezoid<t_vec>>();
			first_bottom->SetLeftPoint(leftpt);
			first_bottom->SetRightPoint(first_trap->GetRightPoint());
			first_bottom->SetTopLine(line);
			first_bottom->SetBottomLine(first_trap->GetBottomLine());


			// last trapezoid
			auto last_trap_node = *intersecting_trapezoids.rbegin();
			auto last_trap = last_trap_node->GetTrapezoid();

			auto last_right = std::make_shared<Trapezoid<t_vec>>();
			last_right->SetLeftPoint(rightpt);
			last_right->SetRightPoint(last_trap->GetRightPoint());
			last_right->SetTopLine(last_trap->GetTopLine());
			last_right->SetBottomLine(last_trap->GetBottomLine());

			auto last_top = std::make_shared<Trapezoid<t_vec>>();
			last_top->SetLeftPoint(last_trap->GetLeftPoint());
			last_top->SetRightPoint(rightpt);
			last_top->SetTopLine(last_trap->GetTopLine());
			last_top->SetBottomLine(line);

			auto last_bottom = std::make_shared<Trapezoid<t_vec>>();
			last_bottom->SetLeftPoint(last_trap->GetLeftPoint());
			last_bottom->SetRightPoint(rightpt);
			last_bottom->SetTopLine(line);
			last_bottom->SetBottomLine(last_trap->GetBottomLine());


			// mid trapezoids
			std::vector<std::shared_ptr<Trapezoid<t_vec>>> mid_tops, mid_bottoms;
			std::vector<std::shared_ptr<TrapezoidNodeTrapezoid<t_vec>>> mid_trap_nodes;

			for(std::size_t isect_idx=1; isect_idx<intersecting_trapezoids.size()-1; ++isect_idx)
			{
				auto mid_trap_node = intersecting_trapezoids[isect_idx];
				auto mid_trap = mid_trap_node->GetTrapezoid();

				auto mid_top = std::make_shared<Trapezoid<t_vec>>();
				mid_top->SetLeftPoint(mid_trap->GetLeftPoint());
				mid_top->SetRightPoint(mid_trap->GetRightPoint());
				mid_top->SetTopLine(mid_trap->GetTopLine());
				mid_top->SetBottomLine(line);

				auto mid_bottom = std::make_shared<Trapezoid<t_vec>>();
				mid_bottom->SetLeftPoint(mid_trap->GetLeftPoint());
				mid_bottom->SetRightPoint(mid_trap->GetRightPoint());
				mid_bottom->SetTopLine(line);
				mid_bottom->SetBottomLine(mid_trap->GetBottomLine());

				if(isect_idx > 1)
				{
					try_unite_trapezoids(&*mid_tops.rbegin(), &mid_top, eps);
					try_unite_trapezoids(&*mid_bottoms.rbegin(), &mid_bottom , eps);
					try_unite_trapezoids(&*mid_tops.rbegin(), &mid_bottom, eps);
					try_unite_trapezoids(&*mid_bottoms.rbegin(), &mid_top, eps);
				}

				if(isect_idx == 1)
				{
					try_unite_trapezoids(&first_top, &mid_top, eps);
					try_unite_trapezoids(&first_bottom, &mid_bottom, eps);
					try_unite_trapezoids(&first_top, &mid_bottom, eps);
					try_unite_trapezoids(&first_bottom, &mid_top, eps);
				}
				else if(isect_idx == intersecting_trapezoids.size()-2)
				{
					try_unite_trapezoids(&last_top, &mid_top, eps);
					try_unite_trapezoids(&last_bottom, &mid_bottom, eps);
					try_unite_trapezoids(&last_top, &mid_bottom, eps);
					try_unite_trapezoids(&last_bottom, &mid_top, eps);
				}

				mid_tops.push_back(mid_top);
				mid_bottoms.push_back(mid_bottom);
				mid_trap_nodes.push_back(mid_trap_node);
			}

			try_unite_trapezoids(&first_top, &last_top, eps);
			try_unite_trapezoids(&first_bottom, &last_bottom, eps);
			try_unite_trapezoids(&first_top, &last_bottom, eps);
			try_unite_trapezoids(&first_bottom, &last_top, eps);


			// first trapezoid
			auto first_left_node = create_trapnode(first_left);
			auto first_top_node = create_trapnode(first_top);
			auto first_bottom_node = create_trapnode(first_bottom);

			auto first_line_node = std::make_shared<TrapezoidNodeLine<t_vec>>(line);
			if(!first_top->IsEmpty(eps))
				first_line_node->SetLeft(first_top_node);
			if(!first_bottom->IsEmpty(eps))
				first_line_node->SetRight(first_bottom_node);

			auto first_leftpt_node = std::make_shared<TrapezoidNodePoint<t_vec>>(leftpt);
			if(!first_left->IsEmpty(eps))
				first_leftpt_node->SetLeft(first_left_node);
			first_leftpt_node->SetRight(first_line_node);

			fit_trapezoid_lines<t_vec>(first_leftpt_node);

			if(root->GetType() == TrapezoidNodeType::TRAPEZOID &&
				first_trap == std::dynamic_pointer_cast<TrapezoidNodeTrapezoid<t_vec>>(root)->GetTrapezoid())
			{
				// replace root node
				root = first_leftpt_node;
			}
			else
			{
				// replace child node pointers
				replace_trapezoid_ptr<t_vec>(root, first_trap_node, first_leftpt_node);
			}


			// mid trapezoids
			for(std::size_t i=0; i<mid_trap_nodes.size(); ++i)
			{
				auto mid_top = mid_tops[i];
				auto mid_bottom = mid_bottoms[i];
				auto mid_trap_node = mid_trap_nodes[i];
				auto mid_trap = mid_trap_node->GetTrapezoid();

				auto mid_top_node = create_trapnode(mid_top);
				auto mid_bottom_node = create_trapnode(mid_bottom);

				auto mid_line_node = std::make_shared<TrapezoidNodeLine<t_vec>>(line);
				if(!mid_top->IsEmpty(eps))
					mid_line_node->SetLeft(mid_top_node);
				if(!mid_bottom->IsEmpty(eps))
					mid_line_node->SetRight(mid_bottom_node);

				fit_trapezoid_lines<t_vec>(mid_line_node);
				//save_trapezoid_svg<t_vec, t_line>(mid_line_node, "mid.svg", nullptr);

				if(root->GetType() == TrapezoidNodeType::TRAPEZOID &&
					mid_trap == std::dynamic_pointer_cast<TrapezoidNodeTrapezoid<t_vec>>(root)->GetTrapezoid())
				{
					// replace root node
					root = mid_line_node;
				}
				else
				{
					// replace child node pointers
					replace_trapezoid_ptr<t_vec>(root, mid_trap_node, mid_line_node);
				}
			}


			// last trapezoid
			auto last_right_node = create_trapnode(last_right);
			auto last_top_node = create_trapnode(last_top);
			auto last_bottom_node = create_trapnode(last_bottom);

			auto last_line_node = std::make_shared<TrapezoidNodeLine<t_vec>>(line);
			if(!last_top->IsEmpty(eps))
				last_line_node->SetLeft(last_top_node);
			if(!last_bottom->IsEmpty(eps))
				last_line_node->SetRight(last_bottom_node);

			auto last_rightpt_node = std::make_shared<TrapezoidNodePoint<t_vec>>(rightpt);
			last_rightpt_node->SetLeft(last_line_node);
			if(!last_right->IsEmpty(eps))
				last_rightpt_node->SetRight(last_right_node);

			fit_trapezoid_lines<t_vec>(last_rightpt_node);

			if(root->GetType() == TrapezoidNodeType::TRAPEZOID &&
				last_trap == std::dynamic_pointer_cast<TrapezoidNodeTrapezoid<t_vec>>(root)->GetTrapezoid())
			{
				// replace root node
				root = last_rightpt_node;
			}
			else
			{
				// replace child node pointers
				replace_trapezoid_ptr<t_vec>(root, last_trap_node, last_rightpt_node);
			}
		}

		// save intermediate steps
		//std::vector<t_line> linetmp{{line}};
		//save_trapezoid_svg<t_vec, t_line>(root, std::string{"step"}+std::to_string(lineidx)+".svg", &linetmp);
	}

	if(shear)
	{
		// TODO: remove tilt
		t_mat shear_inv = m::shear<t_mat>(2, 2, 0, 1, -shear_eps);
		trafo_trapezoid_tree<t_vec, t_mat>(root, shear_inv);
	}

	clean_trapezoid_tree<t_vec>(root);
	return root;
}
// ----------------------------------------------------------------------------


/**
 * split a concave polygon
 * @see algorithm: lecture notes by D. Hegazy, 2015
 */
template<class t_vec, class t_real = typename t_vec::value_type>
requires m::is_vec<t_vec>
std::vector<std::vector<t_vec>> convex_split(
	const std::vector<t_vec>& poly, t_real eps = 1e-6)
{
	bool always_split_on_intersection = false;
	const std::size_t N = poly.size();

	if(N <= 3)
		return {};


	// find concave corner
	std::optional<std::size_t> idx_concave, idx_intersection;

	for(std::size_t idx1=0; idx1<N; ++idx1)
	{
		std::size_t idx2 = (idx1+1) % N;
		std::size_t idx3 = (idx1+2) % N;

		const t_vec& vert1 = poly[idx1];
		const t_vec& vert2 = poly[idx2];
		const t_vec& vert3 = poly[idx3];

		t_real angle = m::pi<t_real>-line_angle<t_vec>(vert1, vert2, vert2, vert3);
		angle = m::mod_pos<t_real>(angle, t_real(2)*m::pi<t_real>);
		//std::cout << "angle: " << angle/m::pi<t_real>*180. << std::endl;

		// corner angle > 180°  =>  concave corner found
		if(angle > m::pi<t_real>)
		{
			idx_concave = idx1;
			break;
		}
	}


	// get intersection of concave edge with contour
	t_vec inters{};

	if(idx_concave)
	{
		std::size_t idx2 = (*idx_concave+1) % N;

		const t_vec& vert1 = poly[*idx_concave];
		const t_vec& vert2 = poly[idx2];
		t_vec dir1 = vert2 - vert1;

		circular_wrapper circularverts(const_cast<std::vector<t_vec>&>(poly));

		auto iterBeg = circularverts.begin() + (*idx_concave + 2);
		auto iterEnd = circularverts.begin() + (*idx_concave + N);

		for(auto iter=iterBeg; iter!=iterEnd; ++iter)
		{
			const t_vec& vert3 = *iter;
			const t_vec& vert4 = *(iter + 1);
			t_vec dir2 = vert4 - vert3;

			// intersect infinite line from concave edge with contour line segment
			auto[pt1, pt2, valid, dist, param1, param2] =
				m::intersect_line_line(vert1, dir1, vert3, dir2, eps);
			inters = pt1;

			if(valid && param2>=0. && param2<1.)
			{
				auto iterInters = (iter+1).GetIter();
				idx_intersection = iterInters - poly.begin();
				break;
			}
		}
	}


	// split polygon along the line [idx_concave+1], [idx_intersection]
	std::vector<std::vector<t_vec>>	split{};
	split.reserve(2);

	if(idx_concave && idx_intersection)
	{
		//std::cout << "split indices: " << *idx_concave << ", " << *idx_intersection << std::endl;
		circular_wrapper circularverts(const_cast<std::vector<t_vec>&>(poly));

		auto iter1 = circularverts.begin() + (*idx_concave);
		auto iter2 = circularverts.begin() + (*idx_intersection);

		std::vector<t_vec> poly1, poly2;

		// sub-polygon 1
		for(auto iter = iter2; true; ++iter)
		{
			poly1.push_back(*iter);;
			if(iter.GetIter() == (iter1+1).GetIter())
				break;
		}

		// sub-polygon 2
		for(auto iter = iter1+1; true; ++iter)
		{
			poly2.push_back(*iter);
			if(iter.GetIter() == (iter2).GetIter())
				break;
		}

		// if there's not enough vertices left, split along the line [idx_concave+1], inters instead
		if(poly1.size() < 3 || poly2.size() < 3 || always_split_on_intersection)
		{
			// insert intersection point before [idx_intersection]
			if(!m::equals<t_vec>(*poly1.begin(), inters, eps))
				poly1.insert(poly1.begin(), inters);

			// exchange [idx_intersection] with intersection point
			if(!m::equals<t_vec>(*std::prev(poly2.end(), 2), inters, eps))
				*std::prev(poly2.end(), 1) = inters;
			else
				poly2.resize(poly2.size()-1);
		}

		// recursively split new polygons
		if(auto subsplit1 = convex_split<t_vec, t_real>(poly1); subsplit1.size())
		{
			for(auto&& newpoly : subsplit1)
				split.emplace_back(std::move(newpoly));
		}
		else
		{
			// poly1 was already convex
			split.emplace_back(std::move(poly1));
		}

		if(auto subsplit2 = convex_split<t_vec, t_real>(poly2); subsplit2.size())
		{
			for(auto&& newpoly : subsplit2)
				split.emplace_back(std::move(newpoly));
		}
		else
		{
			// poly2 was already convex
			split.emplace_back(std::move(poly2));
		}
	}


	return split;
}


}
#endif

/**
 * geometric calculations
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date Oct/Nov-2020
 * @license: see 'LICENSE' file
 *
 * Reference for the algorithms:
 *	- "Algorithmische Geometrie" (2005), ISBN: 978-3540209560 (http://dx.doi.org/10.1007/3-540-27619-X).
 */

#ifndef __GEO2D_ALGOS_H__
#define __GEO2D_ALGOS_H__

#include <vector>
#include <queue>
#include <list>
#include <set>
#include <tuple>
#include <stack>
#include <algorithm>
#include <limits>
#include <random>
#include <iostream>

#include "math_algos.h"
#include "math_conts.h"
#include "geo_conts.h"
#include "helpers.h"

#include <boost/intrusive/bstree.hpp>
#include <boost/intrusive/avltree.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullRidge.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullFacetSet.h>
#include <libqhullcpp/QhullVertexSet.h>


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
	auto [inv_trafo, ok] = m::inv<t_mat>(trafo);
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


template<class t_vec>
std::pair<bool, t_vec>
intersect_lines(const t_vec& pos1a, const t_vec& pos1b,
	const t_vec& pos2a, const t_vec& pos2b, bool only_segments=true)
requires m::is_vec<t_vec>
{
	t_vec dir1 = pos1b - pos1a;
	t_vec dir2 = pos2b - pos2a;

	auto[pt1, pt2, valid, dist, param1, param2] =
		m::intersect_line_line(pos1a, dir1, pos2a, dir2);

	if(!valid)
		return std::make_pair(false, m::create<t_vec>({}));

	if(only_segments && (param1<0. || param1>1. || param2<0. || param2>1.))
		return std::make_pair(false, m::create<t_vec>({}));

	return std::make_pair(true, pt1);
}


template<class t_line /*= std::pair<t_vec, t_vec>*/>
std::tuple<bool, typename t_line::first_type>
intersect_lines(const t_line& line1, const t_line& line2)
requires m::is_vec<typename t_line::first_type> && m::is_vec<typename t_line::second_type>
{
	return intersect_lines<typename t_line::first_type>(
		line1.first, line1.second, line2.first, line2.second, true);
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

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// convex hull algorithms
// Reference: "Algorithmische Geometrie" (2005), ISBN: 978-3540209560, Ch. 4.1, pp. 155f
// ----------------------------------------------------------------------------

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


// tests if the vertex is in the hull
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
// delaunay triangulation
// Reference: "Algorithmische Geometrie" (2005), ISBN: 978-3540209560, Ch. 6, pp. 269f
// ----------------------------------------------------------------------------

/**
 * delaunay triangulation and voronoi vertices
 * @returns [ voronoi vertices, triangles, neighbour triangle indices ]
 */
template<class t_vec>
std::tuple<std::vector<t_vec>, std::vector<std::vector<t_vec>>, std::vector<std::set<std::size_t>>>
calc_delaunay(int dim, const std::vector<t_vec>& verts, bool only_hull)
requires m::is_vec<t_vec>
{
	using namespace m_ops;
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
 * @returns [triangle index, shared index 1, shared index 2, nonshared index]
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
 * see: https://de.wikipedia.org/wiki/Algorithmus_von_Kruskal
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
// Reference: "Algorithmische Geometrie" (2005), ISBN: 978-3540209560, Ch. 4.4, pp. 195f
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
				edges[j].first, edges[j].second, false); ok)
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
// Reference: "Algorithmische Geometrie" (2005), ISBN: 978-3540209560, Ch. 2.3.2, pp. 64f
// ----------------------------------------------------------------------------

template<class t_vec, class t_line = std::pair<t_vec, t_vec>>
std::vector<std::tuple<std::size_t, std::size_t, t_vec>>
intersect_ineff(const std::vector<t_line>& lines)
requires m::is_vec<t_vec>
{
	std::vector<std::tuple<std::size_t, std::size_t, t_vec>> intersections;

	for(std::size_t i=0; i<lines.size(); ++i)
	{
		for(std::size_t j=i+1; j<lines.size(); ++j)
		{
			const t_line& line1 = lines[i];
			const t_line& line2 = lines[j];

			if(auto [intersects, pt] = intersect_lines<t_line>(line1, line2); intersects)
				intersections.emplace_back(std::make_tuple(i, j, pt));
		}
	}

	return intersections;
}


template<class t_hook, class t_vec, class t_line = std::pair<t_vec, t_vec>>
requires m::is_vec<t_vec>
struct IntersTreeLeaf
{
	using t_real = typename t_vec::value_type;

	const t_real *curX{nullptr};
	const std::vector<t_line> *lines{nullptr};
	std::size_t line_idx{0};

	t_hook _h{};

	friend std::ostream& operator<<(std::ostream& ostr, const IntersTreeLeaf<t_hook, t_vec, t_line>& e)
	{
		ostr << std::get<0>((*e.lines)[e.line_idx]) << ", " << std::get<1>((*e.lines)[e.line_idx]);
		return ostr;
	}

	friend bool operator<(const IntersTreeLeaf<t_hook, t_vec, t_line>& e1, const IntersTreeLeaf<t_hook, t_vec, t_line>& e2)
	{
		auto get_line_y = [](const t_line& line, t_real x) -> t_real
		{
			const t_vec& pt1 = std::get<0>(line);
			const t_vec& pt2 = std::get<1>(line);

			t_real slope = (pt2[1]-pt1[1]) / (pt2[0]-pt1[0]);
			return pt1[1] + (x-pt1[0])*slope;
		};

		t_real line1_y = get_line_y((*e1.lines)[e1.line_idx], *e1.curX);
		t_real line2_y = get_line_y((*e2.lines)[e2.line_idx], *e2.curX);

		// compare by y
		return line1_y < line2_y;
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


template<class t_vec, class t_line = std::pair<t_vec, t_vec>, class t_real = typename t_vec::value_type>
std::vector<std::tuple<std::size_t, std::size_t, t_vec>>
intersect_sweep(const std::vector<t_line>& _lines, t_real eps = 1e-6)
requires m::is_vec<t_vec>
{
	using t_mat = m::mat<t_real, std::vector>;

	// look for vertical lines
	bool has_vert_line = 0;
	t_real min_angle_to_y{std::numeric_limits<t_real>::max()};
	std::vector<t_line> lines = _lines;

	for(const t_line& line : lines)
	{
		if(m::equals<t_real>(line.first[0], line.second[0], eps))
		{
			has_vert_line = 1;
		}
		else
		{
			// get angles relative to y axis
			t_real angle_to_y = line_angle<t_vec>(line.first, line.second) + m::pi<t_real>/t_real(2);
			angle_to_y = m::mod_pos<t_real>(angle_to_y, t_real(2)*m::pi<t_real>);
			if(angle_to_y > m::pi<t_real>/t_real(2))
				angle_to_y -= m::pi<t_real>;

			if(std::abs(angle_to_y) < std::abs(min_angle_to_y))
				min_angle_to_y = angle_to_y;
		}
	}

	//if(has_vert_line)
	//	std::cout << "vertical line; next lowest angle: " << min_angle_to_y/m::pi<t_real>*180. << std::endl;

	// rotate all lines
	std::optional<t_mat> rotmat;
	if(has_vert_line)
	{
		rotmat = rotation_2d<t_mat, t_vec>(-min_angle_to_y * t_real(0.5));

		for(t_line& line : lines)
		{
			line.first = *rotmat * line.first;
			line.second = *rotmat * line.second;
		}
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
			using namespace m_ops;

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
				ostr << ", intersection=" << *intersection;
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
		//std::cout << "Event: "; evt.print(std::cout); std::cout << std::endl;

		switch(evt.ty)
		{
			case SweepEventType::LEFT_VERTEX:
			{
				// activate line
				t_leaf *leaf = new t_leaf{.curX=&curX, .lines=&lines, .line_idx=evt.line_idx};
				auto iter = status.insert_equal(*leaf);

				auto iterPrev = (iter == status.begin() ? status.end() : std::prev(iter, 1));
				auto iterNext = std::next(iter, 1);

				// add possible intersection event
				if(iterPrev != iter && iterPrev != status.end())
				{
					//std::cout << "index: " << evt.line_idx << ", prev: " << iterPrev->line_idx << std::endl;

					const t_line& line = lines[evt.line_idx];
					const t_line& linePrev = lines[iterPrev->line_idx];

					if(auto [intersects, pt] = intersect_lines<t_line>(line, linePrev);
					   intersects && !m::equals<t_real>(curX, pt[0], eps))
					{
						SweepEvent evtPrev{.x = pt[0], .ty=SweepEventType::INTERSECTION,
							.lower_idx=iterPrev->line_idx, .upper_idx=evt.line_idx,
							.intersection=pt};
						events.emplace(std::move(evtPrev));
					}
				}

				// add possible intersection event
				if(iterNext != iter && iterNext != status.end())
				{
					//std::cout << "index: " << evt.line_idx << ", next: " << iterNext->line_idx << std::endl;

					const t_line& line = lines[evt.line_idx];
					const t_line& lineNext = lines[iterNext->line_idx];

					if(auto [intersects, pt] = intersect_lines<t_line>(line, lineNext);
					   intersects && !m::equals<t_real>(curX, pt[0], eps))
					{
						SweepEvent evtNext{.x = pt[0], .ty=SweepEventType::INTERSECTION,
							.lower_idx=evt.line_idx, .upper_idx=iterNext->line_idx,
							.intersection=pt};
						events.emplace(std::move(evtNext));
					}
				}

				break;
			}
			case SweepEventType::RIGHT_VERTEX:
			{
				// find current line
				auto iter = std::find_if(status.begin(), status.end(),
					[&evt](const auto& leaf) -> bool
					{ return leaf.line_idx == evt.line_idx; });

				if(iter == status.end())
					continue;

				auto iterPrev = (iter == status.begin() ? status.end() : std::prev(iter, 1));
				auto iterNext = std::next(iter, 1);

				// inactivate current line
				delete &*iter;
				iter = status.erase(iter);


				// add possible intersection event
				if(iterPrev != iterNext && iterPrev != status.end() && iterNext != status.end())
				{
					//std::cout << "prev: " << iterPrev->line_idx << ", next: " << iterNext->line_idx << std::endl;

					const t_line& linePrev = lines[iterPrev->line_idx];
					const t_line& lineNext = lines[iterNext->line_idx];

					if(auto [intersects, pt] = intersect_lines<t_line>(linePrev, lineNext);
					   intersects && !m::equals<t_real>(curX, pt[0], eps))
					{
						SweepEvent evt{.x = pt[0], .ty=SweepEventType::INTERSECTION,
							.lower_idx=iterPrev->line_idx, .upper_idx=iterNext->line_idx,
							.intersection=pt};
						events.emplace(std::move(evt));
					}
				}

				break;
			}
			case SweepEventType::INTERSECTION:
			{
				if(std::find_if(intersections.begin(), intersections.end(),
					[&evt, eps](const auto& inters) -> bool
						{ return m::equals<t_vec>(std::get<2>(inters), *evt.intersection, eps); })
					== intersections.end())
				{
					// report an intersection
					intersections.emplace_back(std::make_tuple(*evt.lower_idx, *evt.upper_idx, *evt.intersection));
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


				std::swap(iterUpper->line_idx, iterLower->line_idx);
				std::swap(iterUpper, iterLower);


				auto iterPrev = (iterUpper == status.begin() ? status.end() : std::prev(iterUpper, 1));
				auto iterNext = std::next(iterLower, 1);


				// add possible intersection event
				if(iterPrev != iterUpper && iterPrev != status.end() && iterUpper != status.end())
				{
					//std::cout << "prev: " << iterPrev->line_idx << ", upper: " << iterUpper->line_idx << std::endl;

					const t_line& linePrev = lines[iterPrev->line_idx];
					const t_line& lineUpper = lines[iterUpper->line_idx];

					if(auto [intersects, pt] = intersect_lines<t_line>(linePrev, lineUpper);
					   intersects && !m::equals<t_real>(curX, pt[0], eps))
					{
						SweepEvent evtPrev{.x = pt[0], .ty=SweepEventType::INTERSECTION,
							.lower_idx=iterPrev->line_idx, .upper_idx=iterUpper->line_idx,
							.intersection=pt};
						events.emplace(std::move(evtPrev));
					}
				}

				// add possible intersection event
				if(iterNext != iterLower && iterNext != status.end() && iterLower != status.end())
				{
					//std::cout << "next: " << iterNext->line_idx << ", lower: " << iterLower->line_idx << std::endl;

					const t_line& lineNext = lines[iterNext->line_idx];
					const t_line& lineLower = lines[iterLower->line_idx];

					if(auto [intersects, pt] = intersect_lines<t_line>(lineNext, lineLower);
					   intersects && !m::equals<t_real>(curX, pt[0], eps))
					{
						SweepEvent evtNext{.x = pt[0], .ty=SweepEventType::INTERSECTION,
							.lower_idx=iterLower->line_idx, .upper_idx=iterNext->line_idx,
							.intersection=pt};
						events.emplace(std::move(evtNext));
					}
				}

				break;
			}
		}
	}


	// rotate intersection points back
	if(rotmat)
	{
		*rotmat = m::trans<t_mat>(*rotmat);

		for(auto& inters : intersections)
			std::get<2>(inters) = *rotmat * std::get<2>(inters);
	}

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
// Reference: "Algorithmische Geometrie" (2005), ISBN: 978-3540209560,
//	Ch. 2.2.2, pp. 53f; Ch. 2.3.1, pp. 57f; Ch. 2.4.1, pp. 93f
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
 * @see:
 *	- http://dx.doi.org/10.1007/3-540-27619-X, ch 2.3.1, p. 57
 *	- https://en.wikipedia.org/wiki/Closest_pair_of_points_problem
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
}
#endif

/**
 * geometric calculations
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 18-Oct-2020
 * @license: see 'LICENSE' file
 *
 * References:
 *	- http://dx.doi.org/10.1007/3-540-27619-X, esp. ch 4.1, ch. 6
 */

#ifndef __GEO2D_H__
#define __GEO2D_H__

#include <vector>
#include <list>
#include <tuple>
#include <algorithm>
#include <limits>

#include "math_algos.h"
#include "math_conts.h"
#include "helpers.h"

#include <boost/intrusive/bstree.hpp>

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullRidge.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>


// ----------------------------------------------------------------------------
// Helper functions
// ----------------------------------------------------------------------------

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
_sort_vertices_by_angle(const std::vector<t_vec>& _verts)
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

// ----------------------------------------------------------------------------


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

		return std::get<0>(_sort_vertices_by_angle<t_vec>(hullverts));
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
	return std::get<0>(_sort_vertices_by_angle<t_vec>(hullLeft));
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
	std::tie(hull, vert_in_hull) = _sort_vertices_by_angle<t_vec>(hull);


	// test if the vertex is already in the hull
	auto is_in_hull = [&hull, &vert_in_hull](const t_vec& newvert) -> std::tuple<bool, std::size_t, std::size_t>
	{
		for(std::size_t hullvertidx1=0; hullvertidx1<hull.size(); ++hullvertidx1)
		{
			std::size_t hullvertidx2 = hullvertidx1+1;
			if(hullvertidx2 >= hull.size())
				hullvertidx2 = 0;

			const t_vec& hullvert1 = hull[hullvertidx1];
			const t_vec& hullvert2 = hull[hullvertidx2];

			// new vertex is between these two points
			if(side_of_line<t_vec>(vert_in_hull, hullvert1, newvert) > 0. &&
				side_of_line<t_vec>(vert_in_hull, hullvert2, newvert) <= 0.)
			{
				// outside hull?
				if(side_of_line<t_vec>(hullvert1, hullvert2, newvert) < 0.)
					return std::make_tuple(false, hullvertidx1, hullvertidx2);
			}
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
	auto is_in_hull = [&vert_in_hull, &hull](const t_vec& newvert) -> std::tuple<bool, std::size_t, std::size_t>
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


/**
 * delaunay triangulation and voronoi vertices
 */
template<class t_vec>
std::tuple<std::vector<t_vec>, std::vector<std::vector<t_vec>>>
calc_delaunay(int dim, const std::vector<t_vec>& verts, bool only_hull)
requires m::is_vec<t_vec>
{
	using namespace m_ops;
	namespace qh = orgQhull;

	using t_real = typename t_vec::value_type;
	using t_real_qhull = coordT;

	std::vector<t_vec> voronoi;		// voronoi vertices
	std::vector<std::vector<t_vec>> triags;	// delaunay triangles

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

		for(auto iterFacet=facets.begin(); iterFacet!=facets.end(); ++iterFacet)
		{
			if(iterFacet->isUpperDelaunay())
				continue;

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

			triags.emplace_back(std::move(thetriag));
		}
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}

	return std::make_tuple(voronoi, triags);
}



/**
 * delaunay triangulation using parabolic trafo
 */
template<class t_vec>
std::tuple<std::vector<t_vec>, std::vector<std::vector<t_vec>>>
calc_delaunay_parabolic(const std::vector<t_vec>& verts)
requires m::is_vec<t_vec>
{
	using namespace m_ops;
	namespace qh = orgQhull;

	using t_real = typename t_vec::value_type;
	using t_real_qhull = coordT;

	const int dim = 2;
	std::vector<t_vec> voronoi;		// voronoi vertices
	std::vector<std::vector<t_vec>> triags;	// delaunay triangles

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

		for(auto iterFacet=facets.begin(); iterFacet!=facets.end(); ++iterFacet)
		{
			if(iterFacet->isUpperDelaunay())
				continue;

			// filter out non-visible part of hull
			qh::QhullHyperplane plane = iterFacet->hyperplane();
			t_vec normal = m::create<t_vec>(dim+1);
			for(int i=0; i<dim+1; ++i)
				normal[i] = t_real{plane[i]};
			// normal pointing upwards?
			if(normal[2] > 0.)
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
			triags.emplace_back(std::move(thetriag));
		}
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}

	return std::make_tuple(voronoi, triags);
}


// ----------------------------------------------------------------------------


#endif

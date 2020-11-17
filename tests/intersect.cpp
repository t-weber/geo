/**
 * line segment intersections
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 11-oct-20
 * @license see 'LICENSE' file
 *
 * References:
 *	- http://dx.doi.org/10.1007/3-540-27619-X, ch 2.3.2, p. 64
 */

#include "../src/geo_algos.h"
using namespace m_ops;


using t_real = double;
using t_vec = m::vec<t_real, std::vector>;
using t_line = std::pair<t_vec, t_vec>;


int main()
{
	std::vector<t_line> lines{{
		std::make_pair(m::create<t_vec>({1., 2.}), m::create<t_vec>({2., 2.})),
		std::make_pair(m::create<t_vec>({1.9, 1.}), m::create<t_vec>({2.1, 3.})),
		std::make_pair(m::create<t_vec>({1.8, 1.1}), m::create<t_vec>({2.2, 3.1})),
		std::make_pair(m::create<t_vec>({0., 0.}), m::create<t_vec>({6., 5.})),
	}};


	{
		auto intersections = intersect_ineff<t_vec, t_line>(lines);
		for(const auto& intersection : intersections)
		{
			std::cout << "Intersection between line " << std::get<0>(intersection)
				<< " and line " << std::get<1>(intersection)
				<< ": " << std::get<2>(intersection)
				<< "." << std::endl;
		}
	}

	std::cout << std::endl;

	{
		auto intersections = intersect_sweep<t_vec, t_line>(lines, 1e-6);
		for(const auto& intersection : intersections)
		{
			std::cout << "Intersection between line " << std::get<0>(intersection)
				<< " and line " << std::get<1>(intersection)
				<< ": " << std::get<2>(intersection)
				<< "." << std::endl;
		}
	}

	return 0;
}

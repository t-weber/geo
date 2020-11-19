/**
 * closest pair
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 4-oct-20
 * @license see 'LICENSE' file
 *
 * References:
 *	- http://dx.doi.org/10.1007/3-540-27619-X, ch 2.3.1, p. 57
 *	- https://en.wikipedia.org/wiki/Closest_pair_of_points_problem
 */

#include "../src/geo_algos.h"


using t_real = double;
using t_vec = m::vec<t_real, std::vector>;
using t_mat = m::mat<t_real, std::vector>;


int main()
{
	using namespace m_ops;

	std::vector<t_vec> points{{
		m::create<t_vec>({1., 0.}),
		m::create<t_vec>({2., 0.5}),
		m::create<t_vec>({3., 7.}),
		m::create<t_vec>({4., 4.}),
		m::create<t_vec>({5., 2.}),
		m::create<t_vec>({6., 3.}),
		m::create<t_vec>({7., 1.}),
		m::create<t_vec>({8., 5.}),
		m::create<t_vec>({9., 5.}),
	}};


	{
		auto [pt1, pt2, dist] = closest_pair_ineff<t_vec>(points);
		if(pt1 && pt2)
		{
			std::cout << "Closest pair (inefficient): point 1: " << *pt1 << ", point 2: " << *pt2
				<< ", dist: " << dist << std::endl;
		}
	}


	{
		auto [pt1, pt2, dist] = closest_pair_sweep<t_vec>(points);
		std::cout << "Closest pair (sweep): point 1: " << pt1 << ", point 2: " << pt2
			<< ", dist: " << dist << std::endl;
	}

	return 0;
}

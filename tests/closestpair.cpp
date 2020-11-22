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


#include <chrono>
#include "../src/geo_algos.h"


using t_real = double;
using t_vec = m::vec<t_real, std::vector>;
using t_mat = m::mat<t_real, std::vector>;


int main()
{
	using namespace m_ops;

	std::size_t num_pts = 1000;
	t_real min = -100;
	t_real max = 100;

	std::vector<t_vec> points;
	points.reserve(num_pts);
	for(std::size_t i=0; i<num_pts; ++i)
		points.emplace_back(m::create<t_vec>({get_rand<t_real>(min, max), get_rand<t_real>(min, max)}));


	{
		auto starttime = std::chrono::steady_clock::now();

		auto [pt1, pt2, dist] = closest_pair_ineff<t_vec>(points);
		if(pt1[0] > pt2[0])
			std::swap(pt1, pt2);

		auto stoptime = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stoptime - starttime);

		if(pt1 && pt2)
		{
			std::cout << "Closest pair (ineff): point 1: " << *pt1 << ", point 2: " << *pt2
				<< ", dist: " << dist << ", time: " << elapsed.count() << " ms" << std::endl;
		}
	}


	{
		auto starttime = std::chrono::steady_clock::now();

		auto [pt1, pt2, dist] = closest_pair_sweep<t_vec>(points);
		//if(pt1[0] > pt2[0])
		//	std::swap(pt1, pt2);

		auto stoptime = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stoptime - starttime);

		std::cout << "Closest pair (sweep): point 1: " << pt1 << ", point 2: " << pt2
			<< ", dist: " << dist << ", time: " << elapsed.count() << " ms" << std::endl;
	}


	{
		auto starttime = std::chrono::steady_clock::now();

		auto [pt1, pt2, dist] = closest_pair_rtree<2, t_vec>(points);
		//if(pt1[0] > pt2[0])
		//	std::swap(pt1, pt2);

		auto stoptime = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stoptime - starttime);

		std::cout << "Closest pair (rtree): point 1: " << pt1 << ", point 2: " << pt2
			<< ", dist: " << dist  << ", time: " << elapsed.count() << " ms" << std::endl;
	}

	return 0;
}

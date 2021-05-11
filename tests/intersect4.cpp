/**
 * line segment intersections
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date may-2021
 * @license see 'LICENSE' file
 */

#include "../src/geo_algos.h"
using namespace m_ops;


using t_real = double;
using t_vec = m::vec<t_real, std::vector>;
using t_line = std::pair<t_vec, t_vec>;


int main()
{
	const t_real min = -1000, max = 1000;

	t_vec pt1a = m::create<t_vec>({g::get_rand(min,max), g::get_rand(min, max)});
	t_vec pt1b = m::create<t_vec>({g::get_rand(min,max), g::get_rand(min, max)});

	t_vec pt2a = m::create<t_vec>({g::get_rand(min,max), g::get_rand(min, max)});
	t_vec pt2b = m::create<t_vec>({g::get_rand(min,max), g::get_rand(min, max)});

	auto [ok, inters] = g::intersect_lines<t_vec>(pt1a, pt1b, pt2a, pt2b, true, 1e-4);
	std::cout << "intersection: " << std::boolalpha << ok << ", at " << inters << "." << std::endl;

	return 0;
}

/**
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date nov-2020
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -Wall -Wextra -Weffc++ -o tst3 tst3.cpp
 */


#include <iostream>
#include <iomanip>
#include <vector>

#include "../src/geo_algos.h"

using t_real = double;
using t_vec = m::vec<t_real, std::vector>;
//using t_mat = m::mat<t_real, std::vector>;


int main()
{
	std::vector<std::pair<t_vec, t_vec>> lines
	{
		// must not cross!
		{ m::create<t_vec>({-5, 5}), m::create<t_vec>({5, 1}) },
		{ m::create<t_vec>({-10, -5}), m::create<t_vec>({10, 0}) },
		{ m::create<t_vec>({-7, -8}), m::create<t_vec>({-1, -3}) },
		{ m::create<t_vec>({6, 7}), m::create<t_vec>({8, 9}) },
		{ m::create<t_vec>({6, 6}), m::create<t_vec>({8, 6}) },
	};

	bool randomise = true;
	bool shear = true;
	auto node = g::create_trapezoid_tree<t_vec>(lines, randomise, shear, 1.);
	std::cout << std::make_pair(node, 0) << std::endl;
	save_trapezoid_svg(node, "tst.svg", &lines);

	return 0;
}

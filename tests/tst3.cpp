/**
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date nov-2020
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -o tst3 tst3.cpp
 */


#include <iostream>
#include <iomanip>
#include <vector>

#include "../src/geo_algos.h"

using t_real = double;
using t_vec = m::vec<t_real, std::vector>;


int main()
{
	std::vector<std::pair<t_vec, t_vec>> lines
	{
		{ m::create<t_vec>({-5, 0}), m::create<t_vec>({5, 0}) },
		{ m::create<t_vec>({-10, 0}), m::create<t_vec>({10, 0}) },
	};

	return 0;
}

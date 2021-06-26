/**
 * split a concave polygon into convex regions
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 26-jun-21
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -Wall -Wextra -o convex convex.cpp
 */

#include "../src/geo_algos.h"
using namespace m_ops;


using t_real = double;
using t_vec = m::vec<t_real, std::vector>;
using t_line = std::pair<t_vec, t_vec>;


int main()
{
	std::vector<t_vec> poly1
	{{
		m::create<t_vec>({0., 0.}),
		m::create<t_vec>({10., 0.}),
		m::create<t_vec>({5., 5.}),
	}};

	std::cout << "poly1" << std::endl;
	auto split1 = g::convex_split<t_vec>(poly1);
	if(split1.size() == 0)
		std::cout << "already convex" << std::endl;
	std::cout << std::endl;


	std::vector<t_vec> poly2
	{{
		m::create<t_vec>({0., 0.}),
		m::create<t_vec>({10., 0.}),
		m::create<t_vec>({5., 5.}),
		m::create<t_vec>({5., 10.}),
	}};

	std::cout << "poly2" << std::endl;
	auto split2 = g::convex_split<t_vec>(poly2);
	if(split2.size() == 0)
		std::cout << "already convex" << std::endl;

	for(std::size_t idx=0; idx<split2.size(); ++idx)
	{
		const auto& polysplit = split2[idx];
		std::cout << "split polygon " << idx << std::endl;

		for(const auto& vec : polysplit)
			std::cout << vec << std::endl;
	}
	std::cout << std::endl;

	return 0;
}

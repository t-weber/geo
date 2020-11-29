/**
 * sorting using convex hull
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date nov-2020
 * @license see 'LICENSE' file
 */


#include <iostream>
#include <iomanip>

#include "geo_algos.h"


using t_real = double;
using t_vec = m::vec<t_real, std::vector>;


int main()
{
	std::size_t num_pts = 100;
	t_real min = -100;
	t_real max = 100;

	std::cout << "Number of random numbers to generate: ";
	std::cin >> num_pts;


	std::vector<t_vec> vecs;

	std::cout << "input sequence: ";
	for(std::size_t i=0; i<num_pts; ++i)
	{
		t_real num = g::get_rand<t_real>(min, max);
		std::cout << num << " ";
		vecs.emplace_back(m::create<t_vec>({ num, num*num }));
	}
	std::cout << "\n" << std::endl;


	auto hull = g::calc_hull_recursive<t_vec>(vecs);


	auto cmp_x = [](const t_vec& vec1, const t_vec& vec2) -> bool
	{
		return vec1[0] < vec2[0];
	};


	// rotate minimum element to front -> O(n)
	auto miniter = std::min_element(hull.begin(), hull.end(), cmp_x);
	std::rotate(hull.begin(), miniter, hull.end());


	std::cout << "sorted output sequence: ";
	for(const auto& vec : hull)
		std::cout << vec[0] << " ";
	std::cout << std::endl;

	bool sorted = std::is_sorted(hull.begin(), hull.end(), cmp_x);
	std::cout << "output is sorted: " << std::boolalpha << sorted << std::endl;

	return 0;
}

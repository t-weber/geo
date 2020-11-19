/**
 * sorting using convex hull
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date nov-2020
 * @license see 'LICENSE' file
 */


#include <iostream>
#include <iomanip>

#include "../src/geo_algos.h"


using t_real = double;
using t_vec = m::vec<t_real, std::vector>;


int main()
{
	std::vector<t_vec> vecs;

	std::cout << "input sequence: ";
	for(std::size_t i=0; i<100; ++i)
	{
		t_real num = get_rand<t_real>(-100, 100);
		std::cout << num << " ";
		vecs.emplace_back(m::create<t_vec>({ num, num*num }));
	}
	std::cout << "\n" << std::endl;


	auto hull = calc_hull_recursive<t_vec>(vecs);


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

/**
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date nov-2020
 * @license see 'LICENSE' file
 */


#include <iostream>
#include <iomanip>
#include <vector>

#include "../src/geo_algos.h"

using t_real = double;
using t_vec = m::vec<t_real, std::vector>;


int main()
{
	t_vec tria1 = m::create<t_vec>({-1, 1});
	t_vec tria2 = m::create<t_vec>({1, 1});
	t_vec tria3 = m::create<t_vec>({0, 3});
	t_vec pt1 = m::create<t_vec>({0, 0});
	t_vec pt2 = m::create<t_vec>({0, 1.5});

	std::cout << "inside 1: " << std::boolalpha << pt_inside_triag(tria1, tria2, tria3, pt1) << std::endl;
	std::cout << "inside 2: " << std::boolalpha << pt_inside_triag(tria1, tria2, tria3, pt2) << std::endl;

	std::cout << std::endl;

	std::cout << "2-norm: " << m::norm<t_vec>(tria1) << std::endl;
	std::cout << "2-norm: " << m::norm<t_vec>(tria1, 2.) << std::endl;
	std::cout << "1-norm: " << m::norm<t_vec>(tria1, 1.) << std::endl;

	return 0;
}

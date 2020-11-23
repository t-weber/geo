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
	t_real x = 2.34;
	for(std::size_t i=0; i<10; ++i)
		std::cout << geo_series(x, i) << std::endl;

	return 0;
}

/**
 * maximum subvector algo test
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 10-apr-20
 * @license see 'LICENSE' file
 *
 * @see https://en.wikipedia.org/wiki/Maximum_subarray_problem
 */

#include <iostream>
#include <iomanip>
#include <array>
#include <type_traits>
#include <cstdint>

#include "../src/geo_algos.h"

//#define DO_TEST


#ifdef DO_TEST

int main()
{
	using t_num = std::int32_t;
	using t_largernum = std::int64_t;
	constexpr std::size_t N = 256;


	for(int i=0; i<1000; ++i)
	{
		std::cout << "\rRun " << i << "   ...   ";
		std::array<t_num, N> deltas;
		std::generate(deltas.begin(), deltas.end(), []() -> t_num
		{
			return get_rand<t_num>();
		});


		auto tup1 = subvec_ineffic<t_largernum>(deltas);
		auto tup2 = subvec_sweep<t_largernum>(deltas);

		if(tup1 == tup2)
		{
			std::cout << "OK" << std::endl;
		}
		else
		{
			std::cerr << "Mismatch: "
				<< "1: " << std::get<0>(tup1) << ", " << std::get<1>(tup1) << ", " << std::get<2>(tup1) << "\n"
				<< "2: " << std::get<0>(tup2) << ", " << std::get<1>(tup2) << ", " << std::get<2>(tup2) << std::endl;
			break;
		}
	}

	return 0;
}


#else


int main()
{
	using t_num = std::int8_t;
	using t_largernum = std::int64_t;
	constexpr std::size_t N = 128;

	std::array<t_num, N> deltas;
	std::generate(deltas.begin(), deltas.end(), []() -> t_num
	{
		return get_rand<t_num>();
	});


	for(std::size_t i=0; i<deltas.size(); ++i)
		std::cout << +deltas[i] << " ";
	std::cout << std::endl;


	{
		auto [start_idx, end_idx, maxval] = subvec_ineffic<t_largernum>(deltas);
		std::cout << "Max. subvec range: [" << start_idx << ", " << end_idx << "[, sum: " << +maxval << std::endl;
	}

	{
		auto [start_idx, end_idx, maxval] = subvec_sweep<t_largernum>(deltas);
		std::cout << "Max. subvec range: [" << start_idx << ", " << end_idx << "[, sum: " << +maxval << std::endl;
	}

	return 0;
}

#endif

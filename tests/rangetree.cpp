/**
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date nov-2020
 * @license see 'LICENSE' file
 */


#include <iostream>
#include <iomanip>
#include <vector>

#include "../src/geo_algos.h"
#include "../src/geo_conts.h"

using t_real = double;
using t_vec = m::vec<t_real, std::vector>;
using t_tree = RangeTree<t_vec>;


int main()
{
	using namespace m_ops;

	t_tree tree;

	tree.insert(m::create<t_vec>( {3., 4.} ));
	tree.insert(m::create<t_vec>( {1., 2.} ));
	tree.insert(m::create<t_vec>( {7., -1.} ));
	tree.insert(m::create<t_vec>( {-5., 10.} ));
	tree.insert(m::create<t_vec>( {-12., 2.} ));
	tree.update();
	std::cout << tree << std::endl;

	std::cout << "\nrange query:\n";
	auto vecs = tree.query_range(m::create<t_vec>( {0., 0.} ), m::create<t_vec>( {10., 10.} ));
	for(const t_vec* vec : vecs)
		std::cout << *vec << "\n";
	std::cout << std::endl;

	return 0;
}

/**
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date nov-2020
 * @license see 'LICENSE' file
 */


#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "../src/geo_algos.h"
#include "../src/geo_conts.h"

using t_real = double;
using t_vec = m::vec<t_real, std::vector>;
using t_tree = g::Treap<t_vec>;


int main()
{
	using namespace m_ops;

	t_tree tree;

	tree.insert(m::create<t_vec>( {1., 6.} ));
	tree.insert(m::create<t_vec>( {2., 8.} ));
	tree.insert(m::create<t_vec>( {3., 7.} ));
	tree.insert(m::create<t_vec>( {4., 9.} ));
	tree.insert(m::create<t_vec>( {5., 10.} ));
	tree.insert(m::create<t_vec>( {6., 1.} ));
	tree.insert(m::create<t_vec>( {7., 2.} ));
	tree.insert(m::create<t_vec>( {8., 5.} ));
	tree.insert(m::create<t_vec>( {9., 4.} ));
	tree.insert(m::create<t_vec>( {10., 3.} ));

	std::ofstream ofstr("treap.graph");
	g::write_graph(ofstr, tree.get_root());

	return 0;
}

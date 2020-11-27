/**
 * geometric container data types
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date Nov-2020
 * @license: see 'LICENSE' file
 *
 * Reference for the algorithms:
 *	- "Algorithmische Geometrie" (2005), ISBN: 978-3540209560 (http://dx.doi.org/10.1007/3-540-27619-X).
 */

#ifndef __GEOCONT_H__
#define __GEOCONT_H__

#include "math_algos.h"
#include "math_conts.h"

#include <boost/intrusive/avltree.hpp>

#include <vector>
#include <iostream>


template<class t_hook, class t_vec>
struct RangeTreeLeaf
{
	using t_real = typename t_vec::value_type;

	const t_vec* vec = nullptr;
	std::size_t idx = 0;
	t_real range = -1.;

	t_hook _h{};


	friend std::ostream& operator<<(std::ostream& ostr, const RangeTreeLeaf<t_hook, t_vec>& e)
	{
		ostr << *e.vec;
		return ostr;
	}

	friend bool operator<(const RangeTreeLeaf<t_hook, t_vec>& e1, const RangeTreeLeaf<t_hook, t_vec>& e2)
	{
		// compare by idx
		return (*e1.vec)[e1.idx] < (*e2.vec)[e2.idx];
	}
};


template<class t_vec>
class RangeTree
{
public:
	using t_node = RangeTreeLeaf<
		boost::intrusive::avl_set_member_hook<
			boost::intrusive::link_mode<
				boost::intrusive::normal_link>>,
		t_vec>;

	using t_tree = boost::intrusive::avltree<
		t_node, boost::intrusive::member_hook<
			t_node, decltype(t_node::_h), &t_node::_h>>;


public:
	RangeTree(std::size_t size)
	{
		m_nodes.reserve(size);
		m_vecs.reserve(size);
	}


	void insert(const t_vec& _vec)
	{
		m_vecs.push_back(_vec);
		const t_vec* vec = &m_vecs[m_vecs.size()-1];

		m_nodes.emplace_back(t_node{.vec=vec, .idx=0});

		m_tree.insert_equal(m_nodes[m_nodes.size()-1]);
	}


private:
	t_tree m_tree{};

	std::vector<t_node> m_nodes{};
	std::vector<t_vec> m_vecs{};
};


#endif

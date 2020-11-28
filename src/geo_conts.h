/**
 * geometric container data types
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date Nov-2020
 * @license: see 'LICENSE' file
 *
 * References:
 *	- "Algorithmische Geometrie" (2005), ISBN: 978-3540209560 (http://dx.doi.org/10.1007/3-540-27619-X).
 *  - https://www.boost.org/doc/libs/1_74_0/doc/html/intrusive/node_algorithms.html
 */

#ifndef __GEOCONT_H__
#define __GEOCONT_H__

/**
 * 1: binary search tree
 * 2: avl tree
 * 3: red-black tree
 */
#ifndef RANGE_TREE_IMPL
	#define RANGE_TREE_IMPL 2
#endif

#include "math_algos.h"
#include "math_conts.h"

#include <iostream>
#include <cstdint>

#if RANGE_TREE_IMPL==1
	#include <boost/intrusive/bstree_algorithms.hpp>
#elif RANGE_TREE_IMPL==2
	#include <boost/intrusive/avltree_algorithms.hpp>
#elif RANGE_TREE_IMPL==3
	#include <boost/intrusive/rbtree_algorithms.hpp>
#endif


template<class t_vec> class RangeTree;


template<class t_vec>
struct RangeTreeNode
{
	using t_real = typename t_vec::value_type;

#if RANGE_TREE_IMPL==1 || RANGE_TREE_IMPL==2
	// bstree or avltree
	using t_balance = std::int64_t;
	t_balance balance = 0;
#elif RANGE_TREE_IMPL==3
	// rbtree
	using t_colour = std::int8_t;
	t_colour colour = 0;
#endif

	RangeTreeNode<t_vec> *parent = nullptr;
	RangeTreeNode<t_vec> *left = nullptr, *right = nullptr;

	// range tree for idx+1
	std::shared_ptr<RangeTree<t_vec>> nextidxtree{};

	std::shared_ptr<const t_vec> vec{};

	std::size_t dim = 0;
	std::size_t idx = 0;
	t_real range[2] = { 0., 0. };


	/**
	 * get all nodes in a linear fashion
	 */
	static void get_vecs(const RangeTreeNode<t_vec>* node, std::vector<std::shared_ptr<const t_vec>>& vecs)
	{
		if(node->left)
			get_vecs(node->left, vecs);

		vecs.push_back(node->vec);

		if(node->right)
			get_vecs(node->right, vecs);
	}


	void print(std::ostream& ostr, std::size_t indent=0) const
	{
		using namespace m_ops;
		const RangeTreeNode<t_vec> *node = this;

		ostr << "ptr: " << (void*)node;
		ostr << ", vec: ";
		if(node->vec)
			ostr << *node->vec;
		else
			ostr << "null";
		ostr << ", idx: " << node->idx;
		ostr << ", range: " << node->range[0] << ".." << node->range[1] << "\n";
		//if(node->nextidxtree)
		//	ostr << "next index tree: " << *node->nextidxtree << "---------------------\n";

		if(node->left || node->right)
		{
			for(std::size_t i=0; i<indent+1; ++i)
				ostr << "  ";
			ostr << "left: ";
			if(node->left)
				node->left->print(ostr, indent+1);
			else
				ostr << "nullptr\n";

			for(std::size_t i=0; i<indent+1; ++i)
				ostr << "  ";
			ostr << "right: ";
			if(node->right)
				node->right->print(ostr, indent+1);
			else
				ostr << "nullptr\n";
		}
	}


	friend std::ostream& operator<<(std::ostream& ostr, const RangeTreeNode<t_vec>& node)
	{
		node.print(ostr);
		return ostr;
	}


	friend bool operator<(const RangeTreeNode<t_vec>& e1, const RangeTreeNode<t_vec>& e2)
	{
		// compare by idx
		return (*e1.vec)[e1.idx] < (*e2.vec)[e2.idx];
	}
};


/**
 * node traits
 * see: https://www.boost.org/doc/libs/1_74_0/doc/html/intrusive/node_algorithms.html
 */
template<class t_vec>
struct RangeTreeNodeTraits
{
	using node = RangeTreeNode<t_vec>;
	using node_ptr = node*;
	using const_node_ptr = const node*;

	// ------------------------------------------------------------------------
#if RANGE_TREE_IMPL==1 || RANGE_TREE_IMPL==2
	// bstree or avltree
	using balance = typename node::t_balance;

	static balance positive() { return 1; }
	static balance negative() { return -1; }
	static balance zero() { return 0; }

	static balance get_balance(const node* thenode)
	{
		if(!thenode) return zero();
		return thenode->balance;
	}

	static void set_balance(node* thenode, balance bal)
	{
		if(!thenode) return;
		thenode->balance = bal;
	}

#elif RANGE_TREE_IMPL==3
	// rbtree
	using color = typename node::t_colour;

	static color red() { return 1; }
	static color black() { return 0; }

	static color get_color(const node* thenode)
	{
		if(!thenode) return black();
		return thenode->colour;
	}

	static void set_color(node* thenode, color col)
	{
		if(!thenode) return;
		thenode->colour = col;
	}
#endif
	// ------------------------------------------------------------------------

	static node* get_parent(const node* thenode)
	{
		if(!thenode) return nullptr;
		return thenode->parent;
	}

	static void set_parent(node* thenode, node* parent)
	{
		if(!thenode) return;
		thenode->parent = parent;
	}

	static node* get_left(const node* thenode)
	{
		if(!thenode) return nullptr;
		return thenode->left;
	}

	static void set_left(node* thenode, node* left)
	{
		if(!thenode) return;
		thenode->left = left;
	}

	static node* get_right(const node* thenode)
	{
		if(!thenode) return nullptr;
		return thenode->right;
	}

	static void set_right(node* thenode, node* right)
	{
		if(!thenode) return;
		thenode->right = right;
	}
};


template<class t_vec>
class RangeTree
{
public:
	using t_node = RangeTreeNode<t_vec>;
	using t_nodetraits = RangeTreeNodeTraits<t_vec>;

#if RANGE_TREE_IMPL==1
	// bstree
	using t_treealgos = boost::intrusive::bstree_algorithms<t_nodetraits>;
#elif RANGE_TREE_IMPL==2
	// avltree
	using t_treealgos = boost::intrusive::avltree_algorithms<t_nodetraits>;
#elif RANGE_TREE_IMPL==3
	// rbtree
	using t_treealgos = boost::intrusive::rbtree_algorithms<t_nodetraits>;
#endif


public:
	RangeTree(std::size_t idx=0) : m_idx{idx}
	{
		t_treealgos::init_header(&m_root);
	}


	~RangeTree()
	{
		free_nodes(root());
	}


	/**
	 * query a rectangular range
	 * TODO
	 */
	std::vector<const t_vec*> query_range(const t_vec& min, const t_vec& max)
	{
		t_node nodemin{.vec=std::make_shared<t_vec>(min), .dim=min.size(), .idx=0};
		t_node nodemax{.vec=std::make_shared<t_vec>(max), .dim=max.size(), .idx=0};

		auto [begin, end] = t_treealgos::bounded_range(&m_root, &nodemin, &nodemax,
			[](const t_node* node1, const t_node* node2) -> bool
			{
				return *node1 < *node2;
			}, true, true);

		std::vector<const t_vec*> vecs;
		for(auto iter=begin; iter!=end; iter=t_treealgos::next_node(iter))
		{
			vecs.push_back(iter->vec.get());
		}

		return vecs;
	}


	void insert(const std::vector<t_vec>& vecs)
	{
		for(const t_vec& vec : vecs)
			insert(vec);
		update();
	}


	void insert(const std::vector<std::shared_ptr<const t_vec>>& vecs)
	{
		for(const std::shared_ptr<const t_vec>& vec : vecs)
			insert(vec);
		update();
	}


	void insert(const t_vec& vec)
	{
		t_node* node = new t_node{.vec=std::make_shared<t_vec>(vec), .dim=vec.size(), .idx=m_idx};
		insert(node);
	}


	void insert(const std::shared_ptr<const t_vec>& vec)
	{
		t_node* node = new t_node{.vec=vec, .dim=vec->size(), .idx=m_idx};
		insert(node);
	}


	void insert(t_node* node)
	{
		t_treealgos::insert_equal(&m_root, root(), node,
			[](const t_node* node1, const t_node* node2) -> bool
			{
				return *node1 < *node2;
			});
	}


	const t_node* root() const
	{
		return t_treealgos::root_node(&m_root);
	}


	t_node* root()
	{
		return t_treealgos::root_node(&m_root);
	}


	void update()
	{
		update(root());
	}


	static void update(t_node* node)
	{
		if(!node) return;

		t_node* left = node->left;
		t_node* right = node->right;

		if(left) update(left);
		if(right) update(right);


		// --------------------------------------------------------------------
		// ranges
		// --------------------------------------------------------------------
		// leaf node
		if(!left && !right && node->vec)
			node->range[0] = node->range[1] = (*node->vec)[node->idx];

		if(left && !right)
		{
			node->range[0] = left->range[0];
			if(node->vec)
				node->range[1] = (*node->vec)[node->idx];
			else
				node->range[1] = left->range[1];
		}

		if(right && !left)
		{
			node->range[1] = right->range[1];
			if(node->vec)
				node->range[0] = (*node->vec)[node->idx];
			else
				node->range[0] = right->range[0];
		}

		if(left && right)
		{
			node->range[0] = left->range[0];
			node->range[1] = right->range[1];
		}
		// --------------------------------------------------------------------


		// --------------------------------------------------------------------
		// subtree for next index
		// --------------------------------------------------------------------
		if(node->idx+1 < node->dim)
		{
			node->nextidxtree = std::make_shared<RangeTree<t_vec>>(node->idx+1);

			std::vector<std::shared_ptr<const t_vec>> vecs;
			RangeTreeNode<t_vec>::get_vecs(node, vecs);

			node->nextidxtree->insert(vecs);
		}
		// --------------------------------------------------------------------
	}


	static void free_nodes(t_node* node)
	{
		if(!node) return;

		free_nodes(node->left);
		free_nodes(node->right);

		delete node;
	}


	friend std::ostream& operator<<(std::ostream& ostr, const RangeTree<t_vec>& tree)
	{
		ostr << *tree.root();
		return ostr;
	}


private:
	t_node m_root{};
	std::size_t m_idx = 0;
};


#endif

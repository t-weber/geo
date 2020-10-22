/**
 * helpers
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date oct-20
 * @license: see 'LICENSE.GPL' file
 */

#ifndef __GEO_HELPERS__
#define __GEO_HELPERS__

#include <iterator>


template<class t_cont>
class circular_iterator : public std::iterator<
	typename std::iterator_traits<typename t_cont::iterator>::iterator_category,
	typename t_cont::value_type,
	typename t_cont::difference_type,
	typename t_cont::pointer,
	typename t_cont::reference>
{
public:
	using t_iter = typename t_cont::iterator;

public:
	circular_iterator(t_cont& cont, t_iter iter, std::size_t round=0)
		: m_cont(cont), m_iter(iter), m_round{round}
	{}

	t_iter GetIter() { return m_iter; }
	std::size_t GetRound() const { return m_round; }

	typename t_cont::reference operator*() { return *m_iter; }


	bool operator==(circular_iterator iter) const
	{ return this->m_iter.operator==(iter.m_iter); }

	bool operator!=(circular_iterator iter) const
	{ return this->m_iter.operator!=(iter.m_iter); }


	circular_iterator& operator++()
	{
		std::advance(m_iter, 1);
		if(m_iter == m_cont.end())
		{
			// wrap around
			m_iter = m_cont.begin();
			++m_round;
		}
		return *this;
	}

	circular_iterator& operator--()
	{
		if(m_iter == m_cont.begin())
		{
			// wrap around
			m_iter = std::prev(m_cont.end(), 1);
			++m_round;
		}
		else
		{
			std::advance(m_iter, -1);
		}
		return *this;
	}


	circular_iterator operator++(int)
	{
		auto curiter = *this;
		this->operator++();
		return curiter;
	}


	circular_iterator operator--(int)
	{
		auto curiter = *this;
		this->operator--();
		return curiter;
	}


	circular_iterator& operator+=(std::size_t num)
	{
		for(std::size_t i=0; i<num; ++i)
			operator++();
		return *this;
	}

	circular_iterator& operator-=(std::size_t num)
	{
		for(std::size_t i=0; i<num; ++i)
			operator--();
		return *this;
	}


	circular_iterator operator+(std::size_t num)
	{
		auto iter = *this;
		iter.operator+=(num);
		return iter;
	}


private:
	t_cont& m_cont;
	t_iter m_iter;

	// how often has the container range been looped?
	std::size_t m_round{0};
};



/**
 * circular access to a container
 */
template<class t_cont>
class circular_wrapper
{
public:
	using iterator = circular_iterator<t_cont>;

public:
	circular_wrapper(t_cont& cont) : m_cont{cont}
	{}


	iterator begin() { return iterator(m_cont, m_cont.begin()); }
	iterator end() { return iterator(m_cont, m_cont.end()); }


	typename t_cont::reference operator[](std::size_t i)
	{ return *(begin() + i); }

	const typename t_cont::reference operator[](std::size_t i) const
	{ return *(begin() + i); }


	iterator erase(iterator begin, iterator end)
	{
		auto idxBeg = begin.GetIter() - m_cont.begin();
		auto idxEnd = end.GetIter() - m_cont.begin();

		if(idxBeg <= idxEnd)
		{
			auto iter = m_cont.erase(std::next(m_cont.begin(),idxBeg), std::next(m_cont.begin(),idxEnd));
			return iterator(m_cont, iter);
		}
		else
		{
			m_cont.erase(std::next(m_cont.begin(),idxBeg), m_cont.end());
			auto iter = m_cont.erase(m_cont.begin(), std::next(m_cont.begin(),idxEnd));
			return iterator(m_cont, iter);
		}
	}

private:
	t_cont& m_cont;
};


#endif

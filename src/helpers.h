/**
 * helpers
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date oct-20
 * @license: see 'LICENSE' file
 */

#ifndef __GEO_HELPERS__
#define __GEO_HELPERS__

//#include <iterator>
#include <tuple>


template<class T, class t_val = typename std::tuple_element<0, T>::type>
constexpr bool is_pair =
	std::is_same_v<T, std::pair<t_val, t_val>> ||
	std::is_same_v<T, std::tuple<t_val, t_val>>;


template<class t_cont>
class circular_iterator /*: public std::iterator<
	typename std::iterator_traits<typename t_cont::iterator>::iterator_category,
	typename t_cont::value_type,
	typename t_cont::difference_type,
	typename t_cont::pointer,
	typename t_cont::reference>*/
{
public:
	using value_type = typename t_cont::value_type;
	using difference_type = typename t_cont::difference_type;
	using pointer = typename t_cont::pointer;
	using reference = typename t_cont::reference;

	using iterator_category = typename std::iterator_traits<typename t_cont::iterator>::iterator_category;
	using t_iter = typename t_cont::/*const_*/iterator;

public:
	circular_iterator(t_cont* cont, t_iter iter, int round=0)
		: m_cont(cont), m_iter(iter), m_round{round}
	{}

	circular_iterator(circular_iterator<t_cont>& iter)
		: m_cont(iter.m_cont), m_iter(iter.m_iter), m_round{iter.m_round}
	{}

	circular_iterator& operator=(circular_iterator<t_cont> iter)
	{
		this->m_cont = iter.m_cont;
		this->m_iter = iter.m_iter;
		this->m_round = iter.m_round;

		return *this;
	}

	t_iter GetIter() { return m_iter; }
	const t_iter GetIter() const { return m_iter; }

	int GetRound() const { return m_round; }
	void SetRound(int round) { m_round = round; }

	typename t_cont::reference operator*() { return *m_iter; }
	typename t_cont::pointer operator->() { return m_iter.operator->(); }


	circular_iterator& operator++()
	{
		std::advance(m_iter, 1);
		if(m_iter == m_cont->end())
		{
			// wrap around
			m_iter = m_cont->begin();
			++m_round;
		}
		return *this;
	}

	circular_iterator& operator--()
	{
		if(m_iter == m_cont->begin())
		{
			// wrap around
			m_iter = std::prev(m_cont->end(), 1);
			--m_round;
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

	circular_iterator operator-(std::size_t num)
	{
		auto iter = *this;
		iter.operator-=(num);
		return iter;
	}


	/*std::ptrdiff_t operator-(circular_iterator iter) const
	{
		return *this - iter;
	}*/


	bool operator==(const circular_iterator& iter) const
	{
		return this->GetRound() == iter.GetRound() &&
			this->GetIter() == iter.GetIter();
	}

	bool operator!=(const circular_iterator& iter) const
	{ return !this->operator==(iter); }


	bool operator<(const circular_iterator& iter) const
	{
		if(this->GetRound() == iter.GetRound())
			return std::distance(iter.GetIter(), this->GetIter()) < 0;
		return this->GetRound() < iter.GetRound();
	}

	bool operator<=(const circular_iterator& iter) const
	{ return this->operator<(iter) || this->operator==(iter); }



private:
	t_cont* m_cont{nullptr};
	t_iter m_iter;

	// how often has the container range been looped?
	int m_round{0};
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
	circular_wrapper(t_cont& cont) : m_cont{&cont}
	{}


	iterator begin() { return iterator(m_cont, m_cont->begin()); }
	iterator end() { return iterator(m_cont, m_cont->end()); }


	typename t_cont::reference operator[](std::size_t i)
	{ return *(begin() + i); }

	const typename t_cont::reference operator[](std::size_t i) const
	{ return *(begin() + i); }


	iterator erase(iterator begin, iterator end)
	{
		std::ptrdiff_t idxBeg = std::distance(m_cont->begin(), begin.GetIter());
		std::ptrdiff_t idxEnd = std::distance(m_cont->begin(), end.GetIter());

		if(idxBeg <= idxEnd && begin.GetRound() == end.GetRound())
		{
			// no wrapping around, simply delete range
			auto iter = m_cont->erase(std::next(m_cont->begin(),idxBeg), std::next(m_cont->begin(), idxEnd));
			return iterator{m_cont, iter};
		}
		else if(idxBeg > idxEnd && end.GetRound() > begin.GetRound())
		{
			// wrapping around, split range into two
			m_cont->erase(std::next(m_cont->begin(),idxBeg), m_cont->end());
			auto iter = m_cont->erase(m_cont->begin(), std::next(m_cont->begin(), idxEnd));
			return iterator{m_cont, iter};
		}

		return end;
	}


private:
	t_cont* m_cont{nullptr};
};


#endif


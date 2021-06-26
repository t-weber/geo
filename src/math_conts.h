/**
 * containers and operators for use with math algorithms
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date jan-2018 - jun-2021
 * @license: see 'LICENSE' file
 */

#ifndef __MATH_CONTS_H__
#define __MATH_CONTS_H__

#include <boost/algorithm/string.hpp>
#include <cassert>
#include <vector>
#include <iostream>
#include <iomanip>
#include "math_concepts.h"


// separator tokens
#define COLSEP ';'
#define ROWSEP '|'


// math operators
namespace m_ops {
// ----------------------------------------------------------------------------
// vector operators
// ----------------------------------------------------------------------------

/**
 * unary +
 */
template<class t_vec>
const t_vec& operator+(const t_vec& vec1)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	return vec1;
}


/**
 * unary -
 */
template<class t_vec>
t_vec operator-(const t_vec& vec1)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	using t_size = decltype(t_vec{}.size());
	t_vec vec(vec1.size());

	for(t_size i=0; i<vec1.size(); ++i)
		vec[i] = -vec1[i];

	return vec;
}


/**
 * binary +
 */
template<class t_vec>
t_vec operator+(const t_vec& vec1, const t_vec& vec2)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	using t_size = decltype(t_vec{}.size());

	if constexpr(m::is_dyn_vec<t_vec>)
		assert((vec1.size() == vec2.size()));
	else
		static_assert(vec1.size() == vec2.size());

	t_vec vec(vec1.size());

	for(t_size i=0; i<vec1.size(); ++i)
		vec[i] = vec1[i] + vec2[i];

	return vec;
}


/**
 * binary -
 */
template<class t_vec>
t_vec operator-(const t_vec& vec1, const t_vec& vec2)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	return vec1 + (-vec2);
}


/**
 * vector * scalar
 */
template<class t_vec>
t_vec operator*(const t_vec& vec1, typename t_vec::value_type d)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	using t_size = decltype(t_vec{}.size());
	t_vec vec(vec1.size());

	for(t_size i=0; i<vec1.size(); ++i)
		vec[i] = vec1[i] * d;

	return vec;
}


/**
 * vector * vector
 */
template<class t_vec>
typename t_vec::value_type operator*(const t_vec& vec1, const t_vec& vec2)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	return inner<t_vec>(vec1, vec2);
}


/**
 * scalar * vector
 */
template<class t_vec>
t_vec operator*(typename t_vec::value_type d, const t_vec& vec)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec> 
	//&& !m::is_basic_mat<typename t_vec::value_type>	// hack!
{
	return vec * d;
}

/**
 * vector / scalar
 */
template<class t_vec>
t_vec operator/(const t_vec& vec, typename t_vec::value_type d)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	using T = typename t_vec::value_type;
	return vec * (T(1)/d);
}


/**
 * vector += vector
 */
template<class t_vec>
t_vec& operator+=(t_vec& vec1, const t_vec& vec2)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	vec1 = vec1 + vec2;
	return vec1;
}

/**
 * vector -= vector
 */
template<class t_vec>
t_vec& operator-=(t_vec& vec1, const t_vec& vec2)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	vec1 = vec1 - vec2;
	return vec1;
}


/**
 * vector *= scalar
 */
template<class t_vec>
t_vec& operator*=(t_vec& vec1, typename t_vec::value_type d)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	vec1 = vec1 * d;
	return vec1;
}

/**
 * vector /= scalar
 */
template<class t_vec>
t_vec& operator/=(t_vec& vec1, typename t_vec::value_type d)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	vec1 = vec1 / d;
	return vec1;
}



/**
 * operator <<
 */
template<class t_vec>
std::ostream& operator<<(std::ostream& ostr, const t_vec& vec)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	using t_size = decltype(t_vec{}.size());
	const t_size N = vec.size();

	for(t_size i=0; i<N; ++i)
	{
		ostr << vec[i];
		if(i < N-1)
			ostr << COLSEP << " ";
	}

	return ostr;
}


/**
 * operator >>
 */
template<class t_vec>
std::istream& operator>>(std::istream& istr, t_vec& vec)
requires m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	vec.clear();

	std::string str;
	std::getline(istr, str);

	std::vector<std::string> vecstr;
	boost::split(vecstr, str, [](auto c)->bool { return c==COLSEP; }, boost::token_compress_on);

	for(auto& tok : vecstr)
	{
		boost::trim(tok);
		std::istringstream istr(tok);

		typename t_vec::value_type c{};
		istr >> c;
		vec.emplace_back(std::move(c));
	}

	return istr;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// matrix operators
// ----------------------------------------------------------------------------

/**
 * unary +
 */
template<class t_mat>
const t_mat& operator+(const t_mat& mat1)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
{
	return mat1;
}


/**
 * unary -
 */
template<class t_mat>
t_mat operator-(const t_mat& mat1)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
{
	using t_size = decltype(t_mat{}.size1());
	t_mat mat(mat1.size1(), mat1.size2());

	for(t_size i=0; i<mat1.size1(); ++i)
		for(t_size j=0; j<mat1.size2(); ++j)
			mat(i,j) = -mat1(i,j);

	return mat;
}


/**
 * binary +
 */
template<class t_mat>
t_mat operator+(const t_mat& mat1, const t_mat& mat2)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
{
	using t_size = decltype(t_mat{}.size1());

	if constexpr(m::is_dyn_mat<t_mat>)
		assert((mat1.size1() == mat2.size1() && mat1.size2() == mat2.size2()));
	else
		static_assert(mat1.size1() == mat2.size1() && mat1.size2() == mat2.size2());

	t_mat mat(mat1.size1(), mat1.size2());

	for(t_size i=0; i<mat1.size1(); ++i)
		for(t_size j=0; j<mat1.size2(); ++j)
			mat(i,j) = mat1(i,j) + mat2(i,j);

	return mat;
}


/**
 * binary -
 */
template<class t_mat>
t_mat operator-(const t_mat& mat1, const t_mat& mat2)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
{
	return mat1 + (-mat2);
}


/**
 * matrix * scalar
 */
template<class t_mat>
t_mat operator*(const t_mat& mat1, typename t_mat::value_type d)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
{
	using t_size = decltype(t_mat{}.size1());
	t_mat mat(mat1.size1(), mat1.size2());

	for(t_size i=0; i<mat1.size1(); ++i)
		for(t_size j=0; j<mat1.size2(); ++j)
			mat(i,j) = mat1(i,j) * d;

	return mat;
}

/**
 * scalar * matrix
 */
template<class t_mat>
t_mat operator*(typename t_mat::value_type d, const t_mat& mat)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
{
	return mat * d;
}

/**
 * matrix / scalar
 */
template<class t_mat>
t_mat operator/(const t_mat& mat, typename t_mat::value_type d)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
{
	using T = typename t_mat::value_type;
	return mat * (T(1)/d);
}


/**
 * matrix-matrix product
 */
template<class t_mat>
t_mat operator*(const t_mat& mat1, const t_mat& mat2)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
{
	using t_size = decltype(t_mat{}.size1());

	if constexpr(m::is_dyn_mat<t_mat>)
		assert((mat1.size2() == mat2.size1()));
	else
		static_assert(mat1.size2() == mat2.size1());

	t_mat matRet(mat1.size1(), mat2.size2());

	for(t_size row=0; row<matRet.size1(); ++row)
	{
		for(t_size col=0; col<matRet.size2(); ++col)
		{
			matRet(row, col) = 0;
			for(t_size i=0; i<mat1.size2(); ++i)
				matRet(row, col) += mat1(row, i) * mat2(i, col);
		}
	}

	return matRet;
}


/**
 * matrix *= scalar
 */
template<class t_mat>
t_mat& operator*=(t_mat& mat1, typename t_mat::value_type d)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
{
	mat1 = mat1 * d;
	return mat1;
}

/**
 * matrix /= scalar
 */
template<class t_mat>
t_mat& operator/=(t_mat& mat1, typename t_mat::value_type d)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
{
	mat1 = mat1 / d;
	return mat1;
}


/**
 * operator <<
 */
template<class t_mat>
std::ostream& operator<<(std::ostream& ostr, const t_mat& mat)
requires m::is_basic_mat<t_mat> //&& m::is_dyn_mat<t_mat>
{
	using t_size = decltype(t_mat{}.size1());

	const t_size ROWS = mat.size1();
	const t_size COLS = mat.size2();

	for(t_size row=0; row<ROWS; ++row)
	{
		for(t_size col=0; col<COLS; ++col)
		{
			ostr << mat(row, col);
			if(col < COLS-1)
				ostr << COLSEP << " ";
		}

		if(row < ROWS-1)
			ostr << ROWSEP << " ";
	}

	return ostr;
}


/**
 * prints matrix in nicely formatted form
 */
template<class t_mat>
std::ostream& niceprint(std::ostream& ostr, const t_mat& mat)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
{
	using t_size = decltype(t_mat{}.size1());

	const t_size ROWS = mat.size1();
	const t_size COLS = mat.size2();

	for(t_size i=0; i<ROWS; ++i)
	{
		ostr << "(";
		for(t_size j=0; j<COLS; ++j)
			ostr << std::setw(ostr.precision()*1.5) << std::right << mat(i,j);
		ostr << ")";

		if(i < ROWS-1)
			ostr << "\n";
	}

	return ostr;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// mixed operators
// ----------------------------------------------------------------------------

/**
 * matrix-vector product
 */
template<class t_mat, class t_vec>
t_vec operator*(const t_mat& mat, const t_vec& vec)
requires m::is_basic_mat<t_mat> && m::is_dyn_mat<t_mat>
	&& m::is_basic_vec<t_vec> && m::is_dyn_vec<t_vec>
{
	using t_size = decltype(t_mat{}.size1());

	if constexpr(m::is_dyn_mat<t_mat>)
		assert((mat.size2() == t_size(vec.size())));
	else
		static_assert(mat.size2() == t_size(vec.size()));


	t_vec vecRet(mat.size1());

	for(t_size row=0; row<mat.size1(); ++row)
	{
		vecRet[row] = typename t_vec::value_type{/*0*/};
		for(t_size col=0; col<mat.size2(); ++col)
		{
			auto elem = mat(row, col) * vec[col];
			vecRet[row] += elem;
		}
	}

	return vecRet;
}
// ----------------------------------------------------------------------------

}


// maths
namespace m {

// ----------------------------------------------------------------------------
// vector and matrix containers
// ----------------------------------------------------------------------------

template<class T=double, template<class...> class t_cont = std::vector>
requires is_basic_vec<t_cont<T>> && is_dyn_vec<t_cont<T>>
class vec : public t_cont<T>
{
public:
	using value_type = T;
	using container_type = t_cont<T>;

	vec() = default;
	vec(std::size_t SIZE) : t_cont<T>(SIZE) {}
	~vec() = default;

	const value_type& operator()(std::size_t i) const { return this->operator[](i); }
	value_type& operator()(std::size_t i) { return this->operator[](i); }

	using t_cont<T>::operator[];

	friend vec operator+(const vec& vec1, const vec& vec2) { return m_ops::operator+(vec1, vec2); }
	friend vec operator-(const vec& vec1, const vec& vec2) { return m_ops::operator-(vec1, vec2); }
	friend const vec& operator+(const vec& vec1) { return m_ops::operator+(vec1); }
	friend vec operator-(const vec& vec1) { return m_ops::operator-(vec1); }

	friend value_type operator*(const vec& vec1, const vec& vec2) { return m_ops::operator*<vec>(vec1, vec2); }
	friend vec operator*(value_type d, const vec& vec1) { return m_ops::operator*(d, vec1); }
	friend vec operator*(const vec& vec1, value_type d) { return m_ops::operator*(vec1, d); }
	friend vec operator/(const vec& vec1, value_type d) { return m_ops::operator/(vec1, d); }

	vec& operator*=(const vec& vec2) { return m_ops::operator*=(*this, vec2); }
	vec& operator+=(const vec& vec2) { return m_ops::operator+=(*this, vec2); }
	vec& operator-=(const vec& vec2) { return m_ops::operator-=(*this, vec2); }
	vec& operator*=(value_type d) { return m_ops::operator*=(*this, d); }
	vec& operator/=(value_type d) { return m_ops::operator/=(*this, d); }

private:
};


template<class T=double, template<class...> class t_cont = std::vector>
requires is_basic_vec<t_cont<T>> && is_dyn_vec<t_cont<T>>
class mat
{
public:
	using value_type = T;
	using container_type = t_cont<T>;

	mat() = default;
	mat(std::size_t ROWS, std::size_t COLS) : m_data(ROWS*COLS), m_rowsize{ROWS}, m_colsize{COLS} {}
	~mat() = default;

	std::size_t size1() const { return m_rowsize; }
	std::size_t size2() const { return m_colsize; }
	const T& operator()(std::size_t row, std::size_t col) const { return m_data[row*m_colsize + col]; }
	T& operator()(std::size_t row, std::size_t col) { return m_data[row*m_colsize + col]; }

	friend mat operator+(const mat& mat1, const mat& mat2) { return m_ops::operator+(mat1, mat2); }
	friend mat operator-(const mat& mat1, const mat& mat2) { return m_ops::operator-(mat1, mat2); }
	friend const mat& operator+(const mat& mat1) { return m_ops::operator+(mat1); }
	friend mat operator-(const mat& mat1) { return m_ops::operator-(mat1); }

	friend mat operator*(const mat& mat1, const mat& mat2) { return m_ops::operator*(mat1, mat2); }
	friend mat operator*(const mat& mat1, value_type d) { return m_ops::operator*(mat1, d); }
	friend mat operator*(value_type d, const mat& mat1) { return m_ops::operator*(d, mat1); }
	friend mat operator/(const mat& mat1, value_type d) { return m_ops::operator/(mat1, d); }

	template<class t_vec> requires is_basic_vec<t_cont<T>> && is_dyn_vec<t_cont<T>>
	friend t_vec operator*(const mat& mat1, const t_vec& vec2) { return m_ops::operator*(mat1, vec2); }

	mat& operator*=(const mat& mat2) { return m_ops::operator*=(*this, mat2); }
	mat& operator+=(const mat& mat2) { return m_ops::operator+=(*this, mat2); }
	mat& operator-=(const mat& mat2) { return m_ops::operator-=(*this, mat2); }
	mat& operator*=(value_type d) { return m_ops::operator*=(*this, d); }
	mat& operator/=(value_type d) { return m_ops::operator/=(*this, d); }

private:
	container_type m_data{};
	std::size_t m_rowsize{}, m_colsize{};
};

// ----------------------------------------------------------------------------

}


#endif

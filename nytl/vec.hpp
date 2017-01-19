// Copyright (c) 2017 nyorain
// Distributed under the Boost Software License, Version 1.0.
// See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt

/// \file A Vector implementation class.

#pragma once

#ifndef NYTL_INCLUDE_VEC
#define NYTL_INCLUDE_VEC

#include <nytl/fwd/vec.hpp> // nytl::Vec typedefs
#include <nytl/vecOps.hpp> // nytl::vec::print

#include <cstddef> // std::size_t
#include <iterator> // std::reverse_iterator
#include <array> // std::array

namespace nytl {

/// \brief Basic Vector template class.
/// Implements the Vector concept from nytl/vecOps and there are all common vector
/// operations defined for it.
/// The default implementation uses stack storage and derives from std::array.
/// \module vec
template<std::size_t D, typename T>
class Vec : public std::array<T, D> {
public:
	using Array = std::array<T, D>;
	using Size = typename Array::size_type;
	using Value = typename Array::value_type;
	using Reference = typename Array::reference;
	using ConstReference = typename Array::const_reference;

	template<Size OD, typename OT> using Rebind = Vec<OD, OT>;
	static constexpr Vec<D, T> create(Size) { return {}; }

	static constexpr auto dim = D;

public:
	template<std::size_t OD, typename OT>
	constexpr explicit operator Vec<OD, OT>() const;
};

// - implementation -
template<std::size_t D, typename T>
template<std::size_t OD, typename OT>
constexpr Vec<D, T>::operator Vec<OD, OT>() const
{
	auto ret = Vec<OD, OT>::create(this->size());
	for(auto i = 0u; i < std::min(ret.size(), this->size()); ++i)
		ret[i] = (*this)[i];
	for(auto i = std::min(ret.size(), this->size()); i < ret.size(); ++i)
		ret[i] = {};
	return ret;
}

// - free operators -
template<std::size_t D, typename T, std::size_t OD, typename OT>
constexpr Vec<D, T>& operator+=(Vec<D, T>& a, const Vec<OD, OT>& b) noexcept
{
	for(std::size_t i = 0; i < std::min(a.size(), b.size()); ++i)
		a[i] += b[i];
	return a;
}

template<std::size_t D, typename T, std::size_t OD, typename OT>
constexpr Vec<D, T>& operator-=(Vec<D, T>& a, const Vec<OD, OT>& b) noexcept
{
	for(std::size_t i = 0; i < std::min(a.size(), b.size()); ++i)
		a[i] -= b[i];
	return a;
}

template<std::size_t D, typename T, typename OT>
constexpr Vec<D, T>& operator*=(Vec<D, T>& vec, OT fac)
{
	for(auto& val : vec)
		val *= fac;
	return vec;
}

template<std::size_t D, typename T, typename OT>
constexpr Vec<D, T>& operator/=(Vec<D, T>& vec, OT fac)
{
	for(auto& val : vec)
		val /= fac;
	return vec;
}

template<std::size_t D1, std::size_t D2, typename T1, typename T2>
constexpr auto operator+(const Vec<D1, T1>& a, const Vec<D2, T2>& b)
{
	auto ret = Vec<std::min(D1, D2), decltype(a[0] + b[0])>::create(std::min(a.size(), b.size()));
	for(auto i = 0u; i < ret.size(); ++i)
		ret[i] = a[i] + b[i];
	return ret;
}

template<std::size_t D1, std::size_t D2, typename T1, typename T2>
constexpr auto operator-(const Vec<D1, T1>& a, const Vec<D2, T2>& b)
{
	auto ret = Vec<std::min(D1, D2), decltype(a[0] - b[0])>::create(std::min(a.size(), b.size()));
	for(auto i = 0u; i < ret.size(); ++i)
		ret[i] = a[i] - b[i];
	return ret;
}

template<std::size_t D, typename T>
constexpr auto operator-(const Vec<D, T>& a)
{
	auto ret = Vec<D, decltype(-a[0])>::create(a.size());
	for(auto i = 0u; i < ret.size(); ++i)
		ret[i] = -a[i];
	return ret;
}

template<std::size_t D, typename F, typename T>
constexpr auto operator*(const F& f, const Vec<D, T>& a)
{
	auto ret = Vec<D, decltype(f * a[0])>::create(a.size());
	for(auto i = 0u; i < ret.size(); ++i)
		ret[i] = f * a[i];
	return ret;
}

template<std::size_t D, typename T1, typename T2>
constexpr auto operator==(const Vec<D, T1>& a, const Vec<D, T2>& b)
{
	for(auto i = 0u; i < a.size(); ++i)
		if(a[i] != b[i]) return false;
	return true;
}

template<std::size_t D, typename T1, typename T2>
constexpr auto operator!=(const Vec<D, T1>& a, const Vec<D, T2>& b)
{
	return !(a == b);
}

template<std::size_t D, typename T>
std::ostream& operator<<(std::ostream& os, const Vec<D, T>& a)
{
	return vec::print(os, a);
}

} // namespace nytl

#endif // header guard

// Copyright (c) 2017 nyorain
// Distributed under the Boost Software License, Version 1.0.
// See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt

/// \file Various vector related operations.

#pragma once

#ifndef NYTL_INCLUDE_VEC_OPS
#define NYTL_INCLUDE_VEC_OPS

#include <nytl/field.hpp> // nytl::FieldTraits
#include <nytl/tmpUtil.hpp> // nytl::templatize
#include <nytl/scalar.hpp> // nytl::accumulate

#include <functional> // std::plus, std::multiplies
#include <stdexcept> // std::invalid_argument
#include <cmath> // std::fma
#include <iosfwd> // std::ostream

namespace nytl::vec {
namespace detail {

/// \brief Helper that asserts that the given vectors have the same dimension.
/// \requires Types 'V1','V2' shall be Vector types.
template<typename V1, typename V2>
struct AssertSameDimensions {
	static constexpr void call(const V1& a, const V2& b)
	{
		if(a.size() != b.size())
			throw std::invalid_argument("nytl::vec: vectors must have same dimension");
	}
};

/// \brief Helper that asserts that a vector of type V1 has dimension Dim.
/// \reuqires Type 'V' shall be a Vector type.
template<unsigned int Dim, typename V>
struct AssertDimension {
	static constexpr void call(const V& a)
	{
		if(a.size() != Dim)
			throw std::invalid_argument("nytl::vec: vector must have specified dimension");
	}
};

/// \brief Asserts that the both given vectors have the same dimension.
/// Will result in a compile time error if possible, otherwise throws
/// std::invalid_argument.
template<typename V1, typename V2>
constexpr void assertSameDimensions(const V1& a, const V2& b)
{
	AssertSameDimensions<V1, V2>::call(a, b);
}

/// \brief Asserts that the givne vector has dimension Dim.
/// Will result in a compile time error if possible, otherwise throws
/// std::invalid_argument.
template<unsigned int Dim, typename V>
constexpr void assertDimension(const V& a)
{
	AssertDimension<Dim, V>::call(a);
}

/// \brief Implements the dot operations for two vectors over the given types.
/// May be specialized by custom vector spaces. The default implements
/// the default dot product definition.
/// \requires Types 'T1','T2' shall represent a field.
/// \requires Types 'V1','V2' shall be Vector types over 'T1','T2'.
template<typename T1, typename T2>
struct Dot {
	template<typename V1, typename V2>
	static constexpr auto call(const V1& a, const V2& b)
	{
		using RetType = decltype(a[0] * b[0] + a[0] * b[0]);
		auto ret = FieldTraits<RetType>::zero;
		for(auto i = 0u; i < a.size(); ++i)
			ret += a[i] * b[i];
		return ret;
	}

	template<typename V1, typename V2>
	static constexpr void check(const V1& a, const V2& b)
	{
		assertDimensions(a, b);
	}
};

/// \brief Computes the angle between two vectors over the given types.
/// May be specialized for custom vectors spaces. The default implements
/// the default angle formula.
/// \requires Types 'T1','T2' shall represent a field.
/// \requires Types 'V1','V2' shall be Vector types over 'T1','T2'.
template<typename T1, typename T2>
struct Angle {
	template<typename V1, typename V2>
	static constexpr auto call(const V1& a, const V2& b)
	{
		using Field = FieldTraits<typename V1::Value>;
		return Field::acos(dot(a, b) / (length(a) * length(b)));
	}

	template<typename V1, typename V2>
	static constexpr void check(const V1& a, const V2& b)
	{
		using Field = FieldTraits<typename V1::Value>;
		assertDimensions(a, b);
		if(length(a) == Field::zero || length(b) == Field::zero)
			throw std::domain_error("nytl::vec::angle: null vector given");
	}
};

/// \brief Computes the cross product of two vectors.
/// May be specialized for custom vectors spaces. The default implements
/// the default cross product formula.
/// \requires Types 'T1','T2' shall represent a field.
/// \requires Types 'V1','V2' shall be Vector types over 'T1','T2'.
template<typename T1, typename T2>
struct Cross {
	template<typename V1, typename V2>
	static constexpr auto call(const V1& a, const V2& b)
	{
		auto ret = V1::template Rebind<3, decltype(a[0] * b[0] - a[0] * b[0])>::create(3);
		ret[0] = (a[1] * b[2]) - (a[2] * b[1]);
		ret[1] = (a[2] * b[0]) - (a[0] * b[2]);
		ret[2] = (a[0] * b[1]) - (a[1] * b[0]);
		return ret;
	}

	template<typename V1, typename V2>
	static constexpr void check(const V1& a, const V2& b)
	{
		assertDimension<3>(a);
		assertDimension<3>(b);
	}
};

} // namespace detail

// Vec operations without argument checking
namespace nocheck {

/// Like dot, but no sanity checks are performed.
template<typename V1, typename V2>
constexpr auto dot(const V1& a, const V2& b)
{
	using DotType = detail::Dot<typename V1::Value, typename V2::Value>;
	return DotType::call(a, b);
}

/// Like angle, but no sanity checks are performed.
/// There might be unexpected results e.g. when both vectors are
/// equal due to rounding errors.
template<typename V1, typename V2>
constexpr auto angle(const V1& a, const V2& b)
{
	using AngleType = detail::Angle<typename V1::Value, typename V2::Value>;
	return AngleType::call(a, b);
}

/// Like cross, but no sanity checks are performed.
template<typename V1, typename V2>
constexpr auto cross(const V1& a, const V2& b)
{
	using CrossType = detail::Cross<typename V1::Value, typename V2::Value>;
	return CrossType::call(a, b);
}

/// Like normalize, but no sanity checks are performed.
template<typename V>
constexpr auto normalize(const V& a)
{
	using Field = FieldTraits<typename V::Value>;
	return (Field::one / length(a)) * a;
}

} // namespace nocheck

/// \brief Sums up all values of the given vector using the + operator.
/// \requires Type 'V' shall be a Vector
/// \module vecOps
template<typename V>
constexpr auto sum(const V& a)
{
	auto zero = FieldTraits<typename V::Value>::zero;
	return accumulate(a.begin(), a.end(), zero, std::plus<>());
}

/// \brief Mutliplies all values of the given vector using the * operator.
/// \requires Type 'V' shall be a non-empty Vector
/// \module vecOps
template<typename V>
constexpr auto multiply(const V& a)
{
	auto one = a[0];
	return accumulate(a.begin() + 1, a.end(), one, std::multiplies<>());
}

/// \brief Calculates the default dot product for the given vectors.
/// Note that this follows the dot definition for real numbers and does
/// not automatically handle the dot definition for complex numbers.
/// \requires Types 'V1' and 'V2' shall be Vectors.
/// \throws std::invalid_argument if the size of the input vectors differs.
/// \module vecOps
template<typename V1, typename V2>
constexpr auto dot(const V1& a, const V2& b)
{
	using DotType = detail::Dot<typename V1::Value, typename V2::Value>;
	DotType::check(a, b);
	return DotType::call(a, b);
}

/// \brief Calculates the angle in radians between two vectors using the dot product.
/// Therefore it will always return the smaller between the both vectors on a
/// plane in which both vectors lay.
/// For two equal vectors, it will return always 0.0.
/// Does only work for real numbers and does not handle complex vectors.
/// \requires Types 'V1', 'V2' shall be Vectors.
/// \throws std::invalid_argument if the size of the input vectors differs.
/// \throws std::domain_error if at lesat one of the given vectors has a length of 0.
/// \module vecOps
template<typename V1, typename V2>
constexpr auto angle(const V1& a, const V2& b)
{
	using AngleType = detail::Angle<typename V1::Value, typename V2::Value>;
	AngleType::check(a, b);
	return AngleType::call(a, b);

	// detail::assertSameDimensions(a, b);
	//
	// auto la = length(a);
	// auto lb = length(b);
	// if(la == Field::zero || lb == Field::zero)
	// 	throw std::domain_error("nytl::vec::angle: both vectors are null");
	//
	// auto res = Field::acos(nocheck::dot(a, b) / (la * lb));
	//
	// // We do this check here to output 0 for angle(a, a).
	// // This might produce nan somtimes due to rounding errors.
	// // res != res is true when res represents nan
	// if(res != res) res = Field::acos(Field::one);
	// return res;
}

/// \brief Calculates the cross product for two 3-dimensional vectors.
/// \requires Types 'V1', 'V2' shall be Vectors over the same 3-dimensional space.
/// \throws std::domain_error if at least on of the input vectors does not have a size of 3.
/// \module vecOps
template<typename V1, typename V2>
constexpr auto cross(const V1& a, const V2& b)
{
	using CrossType = detail::Cross<typename V1::Value, typename V2::Value>;
	CrossType::check(a, b);
	return CrossType::call(a, b);
}

/// \brief Returns the euclidean norm (or length) of the given vector.
/// \requires Type 'V' shall be a Vector.
/// \module vecOps
template<typename V>
constexpr auto length(const V& a)
{
	using Field = FieldTraits<typename V::Value>;
	return Field::sqrt(nocheck::dot(a, a));
}

/// \brief Returns a normalization of the given vector for the euclidean norm.
/// \requires Type 'V' shall be a Vector.
/// \throws std::domain_error if the vector has the length 0.
/// \module vecOps
template<typename V>
constexpr auto normalize(const V& a)
{
	using Field = FieldTraits<typename V::Value>;

	auto la = length(a);
	if(la == Field::zero)
		throw std::domain_error("nytl::vec::normalize: vector has length 0");

	return (Field::one / la) * a;
}

/// \brief Returns the euclidean distance between two vectors.
/// Another way to describe this operation is the length between the
/// difference of the given vectors.
/// \requires Types 'V1','V2' shall be Vector types.
/// \requires The both given vectors shall have the same dimension.
/// Will not check for this but subtract them from each other.
/// \module vecOps
template<typename V1, typename V2>
constexpr auto distance(const V1& a, const V2& b)
{
	return length(a - b);
}

/// \brief Prints the given vector to the given ostream.
/// If this function is used, header <ostream> must be included.
/// This function does not implement operator<< since this operator should only implemented
/// for the Vector implementation types.
/// \requires Type 'V' shall be a Vector
/// \requires There must be an implementation of operator<<(std::ostream&, V::Value).
/// \module vecOps
template<typename V>
std::ostream& print(std::ostream& os, const V& vec)
{
	auto& tos = templatize<V>(os); // we don't want to include ostream
	tos << "(";

	auto it = vec.begin();
	tos << *it;
	while(++it != vec.end())
	tos << ", " << *it;

	tos << ")";
	return os;
}

/// Contains various component-wise operations for Vectors.
namespace cw {

/// \brief Returns a vector holding the component-wise maximum of the given Vectors.
/// \requires Type 'V' shall be a Vector type.
/// \requires The both given vectors shall have the same dimension.
/// \module vecOps
template<typename V>
constexpr auto max(V a, const V& b)
{
	detail::assertSameDimensions(a, b);

	for(auto i = 0u; i < a.size(); ++i)
		if(b[i] > a[i])
			a[i] = b[i];
	return a;
}

/// \brief Returns a vector holding the component-wise maximum of the given Vectors.
/// \requires Type 'V' shall be a Vector type.
/// \requires The both given vectors shall have the same dimension.
/// \module vecOps
template<typename V>
constexpr auto min(V a, const V& b)
{
	detail::assertSameDimensions(a, b);

	for(auto i = 0u; i < a.size(); ++i)
		if(b[i] < a[i])
			a[i] = b[i];
	return a;
}

/// \brief Multiplies the two vectors component wise
/// \requires Types 'V1', 'V2' shall be Vector types over the same space.
/// \module vecOps
template<typename V1, typename V2>
constexpr auto multiply(const V1& a, const V2& b)
{
	detail::assertSameDimensions(a, b);

	auto ret = typename V1::template Rebind<V1::dim, decltype(a[0] * b[0])> {};
	for(auto i = 0u; i < a.size(); ++i)
		ret[i] = a[i] * b[i];
	return ret;
}

/// \brief Component-wise divides the first vector by the second one.
/// Will not perform any zero checks.
/// \requires Types 'V1', 'V2' shall be Vector types over the same space.
/// \module vecOps
template<typename V1, typename V2>
constexpr auto divide(const V1& a, const V2& b)
{
	detail::assertSameDimensions(a, b);

	auto ret = typename V1::template Rebind<V1::dim, decltype(a[0] / b[0])> {};
	for(auto i = 0u; i < a.size(); ++i)
		ret[i] = a[i] / b[i];
	return ret;
}

} // namespace cw
} // namespace nytl::vec

#endif // header guard

// Copyright (c) 2017 nyorain
// Distributed under the Boost Software License, Version 1.0.
// See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt

#pragma once

#ifndef NYTL_INCLUDE_COMPLEX_VEC
#define NYTL_INCLUDE_COMPLEX_VEC

#include <nytl/vecOps.hpp>
#include <complex> // std::complex
#include <cmath> // std::acos

namespace nytl {

/// FullPrecisionField specialization for complex numbers.
template<typename T>
struct FullPrecisionField<std::complex<T>> {
	using type = std::complex<double>;
};

namespace vec::detail {

/// Dot operations specialization for the field of complex numbers.
/// See https://en.wikipedia.org/wiki/Dot_product#Complex_vectors for
/// the reason behind this definition.
template<typename T1, typename T2>
struct Dot<std::complex<T1>, std::complex<T2>> {
	template<typename V1, typename V2>
	static constexpr auto call(const V1& a, const V2& b)
	{
		using RetType = decltype(a[0] * std::conj(b[0]) + a[0] * std::conj(b[0]));
		auto ret = FieldTraits<RetType>::zero;
		for(auto i = 0u; i < a.size(); ++i)
			ret += a[i] * std::conj(b[i]);
		return ret;
	}
};

/// Angle operation specialization for the field of complex numbers.
/// See https://en.wikipedia.org/wiki/Dot_product#Complex_vectors for
/// the reason behind this definition.
template<typename T1, typename T2>
struct Angle<std::complex<T1>, std::complex<T2>> {
	template<typename V1, typename V2>
	static constexpr auto call(const V1& a, const V2& b)
	{
		return std::acos(real(dot(a, b)) / (length(a) * length(b)));
	}
};

} // namespace vec::detail
} // namespace nytl

#endif // header guard

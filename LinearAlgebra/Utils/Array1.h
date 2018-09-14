// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Parallel.h"
#include "../Simd/SimdBlock.h"

#include <vector>
#include <array>

///
/// Custom array with starting index offset.
/// By default, indexing starts at 0.
///
template<typename T>
class alg_array : public std::vector<T>
{
public:
	alg_array(size_t const n = 0) : std::vector<T>(n){}
	alg_array(size_t const n, const T& t) : std::vector<T>(n, t) {}

	ALG_INLINE explicit alg_array(const alg_array<T>& a, const size_t o, const size_t n, const size_t s = 1)
		: std::vector<T>(n)
	{
		auto d = a.data() + o;
		for (size_t i = 0, i0 = 0; i < n; ++i, i0 += s)
		{
			(*this)[i] = d[i0];
		}
	}

#if ALG_RUN_ON_GPU
	alg_array(const alg_array& o)
		: std::vector<T>(o.size())
	{
		for (size_t i = 0; i < o.size(); ++i)
			(*this)[i] = o[i];
	}
	alg_array& operator= (const alg_array& o)
	{
		data();
		this->resize(o.size());
		for (size_t i = 0; i < o.size(); ++i)
			(*this)[i] = o[i];
		return *this;
	}
#endif

	alg_array(const std::initializer_list<T>& data) : std::vector<T>(data) {}

	void fill(const T& t)
	{
		for (auto& a : *this)
		{
			a = t;
		}
	}

#if ALG_RUN_ON_GPU
	mutable alg_gpu_buffer<T> m_buffer;
	alg_gpu_buffer<T> device_data() const
	{
		if (m_buffer == nullptr)
		{
			m_buffer = alg_gpu_buffer<T>(size());
			m_buffer.write_async(std::vector<T>::data(), size());
		}
		return m_buffer;
	}
#endif

#if ALG_RUN_ON_GPU
	const T* data() const
	{
		if (m_buffer() != nullptr)
		{
			m_buffer.read_sync((T*)std::vector<T>::data(), size());
			m_buffer = alg_gpu_buffer<T>();
		}
		return std::vector<T>::data();
	}
	T* data()
	{
		if (m_buffer() != nullptr)
		{
			m_buffer.read_sync(std::vector<T>::data(), size());
			ALG_GPU_CALL(alg_opencl.queue.flush());
			m_buffer = alg_gpu_buffer<T>();
		}
		return std::vector<T>::data();
	}
#endif

#if ALG_RUN_ON_GPU
	ALG_INLINE T& operator[] (const size_t i) { return data()[i]; }
	ALG_INLINE const T& operator[] (const size_t i) const
	{
		return const_cast<const T&>(const_cast<alg_array<T>&>(*this)[i]);
	}
#endif

};	// class alg_array

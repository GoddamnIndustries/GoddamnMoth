// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Parallel.h"

#include <vector>

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once
#include <array>
#include <vector>
#include "Parallel.h"
//#include "../Gpgpu/GpgpuOpenCL.h"

///
/// Custom array with starting index offset.
/// By default, indexing starts at 0.
///
template<typename T>
class alg_array : public std::vector<T>
{
public:
	alg_array(size_t const n = 0) : std::vector<T>(n) {}
	alg_array(size_t const n, const T& t) : std::vector<T>(n, t) {}
	alg_array(const std::initializer_list<T>& data) : std::vector<T>(data) {}

	ALG_INLINE explicit alg_array(const alg_array<T>& a, const size_t o, const size_t n, const size_t s = 1)
		: std::vector<T>(n)
	{
		auto d = a.data() + o;
		for (size_t i = 0, i0 = 0; i < n; ++i, i0 += s)
		{
			(*this)[i] = d[i0];
		}
	}

	ALG_INLINE T& operator[] (const size_t i) { return this->std::vector<T>::operator[](i); }
	ALG_INLINE const T& operator[] (const size_t i) const
	{
		return const_cast<const T&>(const_cast<alg_array<T>&>(*this)[i]);
	}

	void fill(const T& t)
	{
		for (auto& a : *this)
		{
			a = t;
		}
	}

	//alg_gpu_buffer<T> device_data() const;

};	// class alg_array

#if 0
#define ALG_ARRAY_USE_STANDARD_VECTOR 1

#if ALG_ARRAY_USE_STANDARD_VECTOR

template<typename T>
using alg_std_vector = std::vector<T>;

///
/// Dynamic array based on standard vector.
///
template<typename T>
class alg_array : public std::vector<T>
{
public:

	/*///
	/// Initializes an empty array.
	///
	ALG_INLINE alg_array()
		: base<T>()
	{
	}

	///
	/// Copies a slice from other array.
	///
	ALG_INLINE alg_array(const alg_array<T>& other, const size_type offset, const size_type size, const size_type stride = 1)
		: base<T>(size)
	{
		for (size_type i = 0; i < size; ++i)
		{
			(*this)[i] = other[offset + i * stride];
		}
	}*/

	///
	/// Initializes array of a specified size.
	///
	ALG_INLINE explicit alg_array(const size_t size = 0)
		: std::vector<T>(size)
	{
	}

	/*///
	/// Initializes array of a specified size filled with specified values.
	///
	ALG_INLINE explicit alg_array(const size_type size, const_reference fill_with)
		: base<T>(size, fill_with)
	{
	}

	///
	/// Initializes array with a specified data.
	///
	ALG_INLINE explicit alg_array(const std::initializer_list<value_type>& data)
		: base<T>(data)
	{
	}*/

public:

	void fill(const T& t)
	{
		for (auto& a: *this) a = t;
	}

	ALG_INLINE T& operator[] (const size_t i) { return std::vector<T>::operator[](i); }
	ALG_INLINE const T& operator[] (const size_t i) const
	{
		return const_cast<const T&>(const_cast<alg_array<T>&>(*this)[i]);
	}

};  // class alg_array

#else

///
/// Random-access iterator for indexed objects.
///
template<typename U>
class alg_index_iterator
{
	using value_type = typename U::value_type;
	using size_type = typename U::size_type;
	using difference_type = typename U::difference_type;
	using reference = algstd::conditional_t<algstd::is_const_v<U>, typename U::const_reference, typename U::reference>;
	using pointer = algstd::conditional_t<algstd::is_const_v<U>, typename U::const_pointer, typename U::pointer>;

private:
	U* m_array;
	size_type  m_index;

public:
	alg_index_iterator(U* const u = nullptr, const size_t i = 0) : m_array(u), m_index(i) {}
};  // class iterator

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++                                  Arrays.                                 +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

///
/// Dynamic array with reference counting and slicing support.
///
template<typename T>
class alg_array
{
	using this_type = alg_array<T>;
public:
	using value_type = T;
	using allocator_type = void;
	using reference = T&;
	using const_reference = const T&;
	using pointer = T*;
	using const_pointer = const T*;
	using iterator = alg_index_iterator<this_type>;
	using const_iterator = alg_index_iterator<const this_type>;
	using difference_type = ptrdiff_t;
	using size_type = size_t;

private:
	using data_counter = std::atomic_size_t;
	using data_pointer = struct data_type {
		data_counter counter;
		value_type data[1];
	}*;

private:
	data_pointer m_outer_data;
	size_type m_outer_size;
	size_type m_inner_size;
	size_type m_inner_offset;
	size_type m_inner_stride;

public:

	// ------------------------------------------------------------------------------------ //
	// -----                    Construction and reference counting.                  ----- //
	// ------------------------------------------------------------------------------------ //

	///
	/// Initializes an empty array.
	///
	ALG_INLINE alg_array()
		: m_outer_data(nullptr), m_outer_size(0)
		, m_inner_size(0), m_inner_offset(0), m_inner_stride(1)
	{
	}

	///
	/// Moves other array.
	///
    ALG_INLINE alg_array(alg_array<T>&& other)
		: m_outer_data(other.m_outer_data), m_outer_size(other.m_outer_size)
        , m_inner_size(other.m_inner_size), m_inner_offset(other.m_inner_offset), m_inner_stride(other.m_inner_stride)
	{
		other.m_outer_data = nullptr;
		other.m_outer_size = other.m_inner_size = other.m_inner_offset = 0;
		other.m_inner_stride = 1;
	}

	///
	/// Copies other array.
	/// @note This operation is inexpensive due to arrays are reference-counted.
	///
	ALG_INLINE alg_array(const alg_array<T>& other)
		: m_outer_data(other.m_outer_data), m_outer_size(other.m_outer_size)
		, m_inner_size(other.m_inner_size), m_inner_offset(other.m_inner_offset), m_inner_stride(other.m_inner_stride)
	{
		add_ref();
	}

	///
	/// Copies a slice from other array.
	/// @note This operation is inexpensive due to arrays are reference-counted.
	///
	ALG_INLINE alg_array(const alg_array<T>& other, const size_type offset, const size_type size, const size_type stride = 1)
		: m_outer_data(other.m_outer_data), m_outer_size(other.m_outer_size)
		, m_inner_size(size), m_inner_offset(other.m_inner_offset + offset), m_inner_stride(other.m_inner_stride * stride)
	{
		add_ref();
	}

	///
	/// Initializes array of a specified size.
	///
	ALG_INLINE explicit alg_array(const size_type size)
		: m_outer_data(alloc_construct(size)), m_outer_size(size)
		, m_inner_size(size), m_inner_offset(0), m_inner_stride(1)
	{
	}

	///
	/// Initializes array of a specified size filled with specified values.
	///
	ALG_INLINE explicit alg_array(const size_type size, const_reference fill_with)
		: m_outer_data(alloc_construct(size, fill_with)), m_outer_size(size)
		, m_inner_size(size), m_inner_offset(0), m_inner_stride(1)
	{
	}

	///
	/// Initializes array with a specified data.
	///
	ALG_INLINE explicit alg_array(const std::initializer_list<value_type>& data)
		: m_outer_data(alloc_construct(data.size(), data.begin())), m_outer_size(data.size())
		, m_inner_size(data.size()), m_inner_offset(0), m_inner_stride(1)
	{
	}

	///
	/// Deletes this reference to the array data.
	/// If current instance is last, deallocates the data.
	///
	ALG_INLINE ~alg_array()
	{
		release();
	}

public:

	///
	/// Move-assigns other array to this.
	///
	ALG_INLINE alg_array<T>& operator= (alg_array<T>&& other)
	{
		if (this != &other)
		{
			release();
			m_outer_data = other.m_outer_data;
			m_outer_size = other.m_inner_size;
			m_inner_size = other.m_inner_size;
			m_inner_offset = other.m_inner_offset;
			m_inner_stride = other.m_inner_stride;

			other.m_outer_data = nullptr;
			other.m_outer_size = other.m_inner_size = other.m_inner_offset = 0;
			other.m_inner_stride = 1;
		}
		return *this;
	}

	///
	/// Assigns other array to this.
	///
	ALG_INLINE alg_array<T>& operator= (const alg_array<T>& other)
	{
		if (this != &other)
		{
			release();
			m_outer_data = other.m_outer_data;
			m_outer_size = other.m_inner_size;
			m_inner_size = other.m_inner_size;
			m_inner_offset = other.m_inner_offset;
			m_inner_stride = other.m_inner_stride;
			add_ref();
		}
		return *this;
	}

private:

	///
	/// Allocates data structure.
	///
	ALG_INLINE static data_pointer alloc(const size_type size)
	{
		ALG_ASSERT(size > 0);

		// Allocating memory:
		const data_pointer pointer = static_cast<data_pointer>(::operator new(
			sizeof(data_type) + (size - 1) * sizeof(value_type)));
		if (pointer == nullptr)
		{
			throw std::bad_alloc();
		}
		return pointer;
	}

	///
	/// Allocates and constructs data structure with reference counter set to 1.
	///
	ALG_INLINE static data_pointer alloc_construct(const size_type size)
	{
		if (size != 0)
		{
			const data_pointer pointer = alloc(size);
			new (&pointer->counter) data_counter(1);
			/*for (size_type i = 0; i < size; ++i)
			{
				new (pointer->data + i) value_type();
			}*/
			return pointer;
		}
		return nullptr;
	}

	///
	/// Allocates and constructs data structure filling it with specified value and setting reference counter to 1.
	///
	ALG_INLINE static data_pointer alloc_construct(const size_type size, const_reference fill_with)
	{
		if (size != 0)
		{
			const data_pointer pointer = alloc(size);
			new (&pointer->counter) data_counter(1);
			for (size_type i = 0; i < size; ++i)
			{
				new (pointer->data + i) value_type(fill_with);
			}
			return pointer;
		}
		return nullptr;
	}

	///
	/// Allocates and constructs data structure copying data from iterators and setting reference counter to 1.
	///
	template<typename TIterator>
	ALG_INLINE static data_pointer alloc_construct(const size_type size, const TIterator iter)
	{
		if (size != 0)
		{
			const data_pointer pointer = alloc(size);
			new (&pointer->counter) data_counter(1);
			for (size_type i = 0; i < size; ++i)
			{
				new (pointer->data + i) value_type(*(iter + i));
			}
			return pointer;
		}
		return nullptr;
	}

	///
	/// De-allocates data destructs structure.
	///
	ALG_INLINE static void dealloc_destruct(const data_pointer pointer, const size_type size)
	{
		if (size == 0)
		{
			// For zero-sized data pointers should be zero.
			ALG_ASSERT(pointer == nullptr);
		}
		else
		{
			ALG_ASSERT(pointer != nullptr);

			// Deinitializing and deallocating data.
			for (size_type i = 0; i < size; ++i)
			{
				pointer->data[i].~value_type();
			}
			pointer->counter.~data_counter();
			::operator delete(pointer);
		}
	}

private:

	///
	/// Increments a value of the reference counter.
	///
	void add_ref()
	{
		if (m_outer_data != nullptr)
		{
			++m_outer_data->counter;
		}
	}

	///
	/// Decrements value of the reference counter and destructs data if counter reaches zero.
	///
	void release()
	{
		if (m_outer_data != nullptr)
		{
			if (--m_outer_data->counter == 0)
			{
				dealloc_destruct(m_outer_data, m_outer_size);
				m_outer_data = nullptr;
				m_outer_size = m_inner_size = 0;
			}
		}
	}

public:

	// ------------------------------------------------------------------------------------ //
	// -----                       Internal data structure access.                    ----- //
	// ------------------------------------------------------------------------------------ //

	/// Return iterator referring to the first element of the array.
	/// @{
	ALG_INLINE iterator begin() { return iterator(this); }
	ALG_INLINE const_iterator begin() const { return const_iterator(this); }
	ALG_INLINE const_iterator cbegin() const { return begin(); }
	/// @}

	///
	/// Returns an iterator referring to the past-the-end element of the array.
	/// @{
	ALG_INLINE iterator end() { return iterator(this, m_inner_size); }
	ALG_INLINE const_iterator end() const { return const_iterator(this, m_inner_size); }
	ALG_INLINE const_iterator cend() const { return end(); }
	/// @}

public:

	///
	/// Returns pointer to the data stored in this array slice.
	/// @note Data may not be contiguously stored at the pointer.
	/// @{
	ALG_INLINE pointer data()
	{
		return m_outer_data != nullptr ? m_outer_data->data + m_inner_offset : nullptr;
	}
	ALG_INLINE const_pointer data() const
	{
		return const_cast<const_pointer>(const_cast<this_type*>(this)->data());
	}
	/// @}

	///
	/// Returns stride of data of this array slice.
	///
	ALG_INLINE size_type stride() const
	{
		return m_inner_stride;
	}

	///
	/// Returns size of this array slice.
	///
	ALG_INLINE size_type size() const
	{
		return m_inner_size;
	}

public:

	///
	/// Returns reference to the element of this array slice.
	/// @{
	ALG_INLINE T& operator[] (const size_t i)
	{
		//ALG_ASSERT(i < m_inner_size);
		//return this->data()[i * m_inner_stride];
		return this->m_outer_data->data[i];
	}
	ALG_INLINE const T& operator[] (const size_t i) const
	{
		return const_cast<const T&>(const_cast<this_type&>(*this)[i]);
	}
	/// @}

private:

	// ------------------------------------------------------------------------------------ //
	// -----                    Internal data structure modification.                 ----- //
	// ------------------------------------------------------------------------------------ //

	/// Makes a resized copy of this array slice.
	ALG_INLINE alg_array<T> resize_copy(const size_t size, const_reference fill_rest_with = value_type()) const
	{
		alg_array<T> copy(size, fill_rest_with);
		for (size_type i = 0; i < std::min(m_inner_size, size); ++i)
		{
			copy[i] = (*this)[i];
		}
		return copy;
	}

	/// Makes a resized copy of this array slice.
	ALG_INLINE alg_array<T> resize_copy_deep(const size_t size, const_reference fill_rest_with = value_type()) const
	{
		alg_array<T> copy(size, fill_rest_with);
		for (size_type i = 0; i < std::min(m_inner_size, size); ++i)
		{
			copy[i] = T((*this)[i].copy());
		}
		return copy;
	}

public:

	///
	/// Makes a resized copy of this array slice and assigns to itself.
	///
	ALG_INLINE void resize(const size_t size, const_reference fill_rest_with = value_type())
	{
		*this = std::move(resize_copy(size, fill_rest_with));
	}

	///
	/// Makes a resized copy of this array slice and assigns to itself.
	///
	ALG_INLINE void resize_deep(const size_t size, const_reference fill_rest_with = value_type())
	{
		*this = std::move(resize_copy_deep(size, fill_rest_with));
	}

	///
	/// Makes a resized copy of this array slice and assigns to itself.
	///
	ALG_INLINE alg_array<T> copy() const
	{
		return resize_copy(m_inner_size);
	}

	///
	/// Makes a resized copy of this array slice and assigns to itself.
	/// Elements are also copied.
	///
	ALG_INLINE alg_array<T> deep_copy() const
	{
		return resize_copy_deep(m_inner_size);
	}

	void fill(const T& t)
	{
		for (size_type i = 0; i < m_inner_size; ++i)
		{
			(*this)[i] = t;
		}
	}

	void fill_deep(const T& t)
	{
		for (size_type i = 0; i < m_inner_size; ++i)
		{
			(*this)[i] = T(t.copy());
		}
	}

};	// class alg_array

#endif
#endif
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <new>
#include <array>
#include <vector>
#include <atomic>

#include "Parallel.h"

enum alg_data_location
{
	alg_data_location_host,
	alg_data_location_device,
};	// enum alg_array_location
template<typename T, alg_data_location TLoc>
class alg_array_impl;

///
/// Dynamic array with reference counting and slicing support.
/// Data is forced to be stored in system memory.
///
template<typename T>
using alg_host_array = alg_array_impl<T, alg_data_location_host>;

///
/// Dynamic array with reference counting and slicing support. 
/// Data is forced to be stored on the GPU device.
///
template<typename T>
using alg_device_array = alg_array_impl<T, alg_data_location_device>;

///
/// Dynamic array with reference counting and slicing support.
/// Scalar types are stored on the GPU device, others in system memory.
///
template<typename T>
using alg_array = alg_host_array<T>;//algstd::conditional_t<algstd::is_floating_point_v<T>, alg_device_array<T>, alg_host_array<T>>;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++                           Dynamic Host Arrays.                           +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

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

///
/// Dynamic array with reference counting and slicing support.
///
template<typename T>
class alg_array_impl<T, alg_data_location_host>
{
	using this_type = alg_array_impl<T, alg_data_location_host>;
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
	ALG_INLINE alg_array_impl()
		: m_outer_data(nullptr), m_outer_size(0)
		, m_inner_size(0), m_inner_offset(0), m_inner_stride(1)
	{
	}

	///
	/// Moves other array.
	///
    ALG_INLINE alg_array_impl(this_type&& other) noexcept
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
	ALG_INLINE alg_array_impl(const this_type& other)
		: m_outer_data(other.m_outer_data), m_outer_size(other.m_outer_size)
		, m_inner_size(other.m_inner_size), m_inner_offset(other.m_inner_offset), m_inner_stride(other.m_inner_stride)
	{
		add_ref();
	}

	///
	/// Copies a slice from other array.
	/// @note This operation is inexpensive due to arrays are reference-counted.
	///
	ALG_INLINE alg_array_impl(const this_type& other, const size_type offset, const size_type size, const size_type stride = 1)
		: m_outer_data(other.m_outer_data), m_outer_size(other.m_outer_size)
		, m_inner_size(size), m_inner_offset(other.m_inner_offset + offset), m_inner_stride(other.m_inner_stride * stride)
	{
		add_ref();
	}

	///
	/// Initializes array of a specified size.
	///
	ALG_INLINE explicit alg_array_impl(const size_type size)
		: m_outer_data(alloc_construct(size)), m_outer_size(size)
		, m_inner_size(size), m_inner_offset(0), m_inner_stride(1)
	{
	}

	///
	/// Initializes array of a specified size filled with specified values.
	///
	ALG_INLINE explicit alg_array_impl(const size_type size, const_reference fill_with)
		: m_outer_data(alloc_construct(size, fill_with)), m_outer_size(size)
		, m_inner_size(size), m_inner_offset(0), m_inner_stride(1)
	{
	}

	///
	/// Initializes array with a specified data.
	///
	ALG_INLINE explicit alg_array_impl(const std::initializer_list<value_type>& data)
		: m_outer_data(alloc_construct(data.size(), data.begin())), m_outer_size(data.size())
		, m_inner_size(data.size()), m_inner_offset(0), m_inner_stride(1)
	{
	}

	///
	/// Deletes this reference to the array data.
	/// If current instance is last, deallocates the data.
	///
	ALG_INLINE ~alg_array_impl()
	{
		release();
	}

public:

	///
	/// Move-assigns other array to this.
	///
	ALG_INLINE this_type& operator= (this_type&& other) noexcept
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
	ALG_INLINE this_type& operator= (const this_type& other)
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
	ALG_INLINE this_type resize_copy(const size_t size, const_reference fill_rest_with = value_type()) const
	{
		this_type copy(size, fill_rest_with);
		for (size_type i = 0; i < std::min(m_inner_size, size); ++i)
		{
			copy[i] = (*this)[i];
		}
		return copy;
	}

	/// Makes a resized copy of this array slice.
	ALG_INLINE this_type resize_copy_deep(const size_t size, const_reference fill_rest_with = value_type()) const
	{
		this_type copy(size, fill_rest_with);
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
	ALG_INLINE this_type copy() const
	{
		return resize_copy(m_inner_size);
	}

	///
	/// Makes a resized copy of this array slice and assigns to itself.
	/// Elements are also copied.
	///
	ALG_INLINE this_type deep_copy() const
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

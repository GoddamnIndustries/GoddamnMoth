// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// A very simple linear algebra library and linear system solving library.
// Copyright (C) 2017 Butakov Oleg.
// All rights reserved.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <vector>
#include "cl.hpp"
#include "../Scalar.h"

#define ALG_IS_GPGPU_AVAILABLE_FOR_T(T) (algstd::is_floating_point_v<T>)
#define ALG_GPU_CALL(...) ALG_ASSERT((__VA_ARGS__) == 0)

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// +++++                       OpenCL-accelerated operations.                     +++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

///
///
///
struct alg_opencl_t
{
	cl::Context context;
	cl::Device default_device;
	cl::CommandQueue queue;

	alg_opencl_t()
	{
		std::vector<cl::Platform> all_platforms;
		ALG_GPU_CALL(cl::Platform::get(&all_platforms));
		if (all_platforms.size() == 0) {
			std::cerr << " No platforms found. Check OpenCL installation!\n";
			ALG_ASSERT(0);
		}

		cl::Platform default_platform = all_platforms.front();
		std::cerr << "Using platform: " << default_platform.getInfo<CL_PLATFORM_NAME>() << "\n";

		std::vector<cl::Device> all_devices;
		ALG_GPU_CALL(default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices));
		if (all_devices.size() == 0) {
			std::cerr << " No devices found. Check OpenCL installation!\n";
			ALG_ASSERT(0);
		}
		default_device = all_devices.back();
		std::cerr << "Using device: " << default_device.getInfo<CL_DEVICE_NAME>() << "\n";
		std::cerr << "Units: " << default_device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << "\n";

		context = cl::Context(all_devices);
		queue = cl::CommandQueue(context, default_device);
	}

} static alg_opencl;

// ------------------------------------------------------------------------------------ //
// -----                               GPU buffers.                               ----- //
// ------------------------------------------------------------------------------------ //

#define ALG_GPU_REDUCE_BLOCK_SIZE 4
#define ALG_GPU_LOCAL_BLOCK_SIZE 128

ALG_INLINE size_t alg_divide_ceil(const size_t numerator, const size_t denominator)
{
	return (size_t)std::ceil((float)numerator / denominator);
}
ALG_INLINE size_t alg_divide_ceil_mul(const size_t numerator, const size_t denominator)
{
	return alg_divide_ceil(numerator, denominator) * denominator;
}

using alg_gpu_event = cl::Event;

///
/// Reference-counted GPU buffer.
///
template<typename T>
struct alg_gpu_buffer : public cl::Buffer
{
public:
	ALG_INLINE explicit alg_gpu_buffer() {}
	ALG_INLINE explicit alg_gpu_buffer(const cl::Buffer& buffer) : cl::Buffer(buffer) {}
	ALG_INLINE explicit alg_gpu_buffer(const size_t size, const T* const data = nullptr)
    {
	    cl_int result;
		cl::Buffer::operator=(cl::Buffer(alg_opencl.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, size * sizeof(T),
			const_cast<T*>(data), &result));
	    ALG_GPU_CALL(result);
    }

public:
	ALG_INLINE bool operator== (std::nullptr_t) const { return (*this)() == nullptr; }
	ALG_INLINE bool operator!= (std::nullptr_t) const { return (*this)() != nullptr; }

public:
	ALG_INLINE void write_async(const T* const data, const size_t size)
	{
		ALG_ASSERT(data != nullptr);
		ALG_GPU_CALL(alg_opencl.queue.enqueueWriteBuffer(*this, CL_TRUE, 0, size * sizeof(T), data));
	}
	ALG_INLINE void write_sync(const T* const data, const size_t size)
	{
		write_async(data, size);
		ALG_GPU_CALL(alg_opencl.queue.flush());
	}

	ALG_INLINE void read_async(T* const data, const size_t size) const
	{
		ALG_ASSERT(data != nullptr);
		ALG_GPU_CALL(alg_opencl.queue.enqueueReadBuffer(*this, CL_TRUE, 0, size * sizeof(T), data));
	}
	ALG_INLINE void read_sync(T* const data, const size_t size) const
	{
		read_async(data, size);
		ALG_GPU_CALL(alg_opencl.queue.flush());
	}
};  // struct alg_gpu_buffer<T>

///
/// Reference-counted sized GPU buffer.
///
template<typename T>
struct alg_gpu_sized_buffer : alg_gpu_buffer<T>
{
	size_t size;

public:
	ALG_INLINE alg_gpu_sized_buffer(const cl::Buffer& buffer, const size_t size)
		: alg_gpu_buffer<T>(buffer), size(size) {}

public:
	ALG_INLINE void write_async(const T* const data)
	{
		alg_gpu_buffer<T>::write_async(data, size);
	}
	ALG_INLINE void write_sync(const T* const data)
	{
		alg_gpu_buffer<T>::write_sync(data, size);
	}

	ALG_INLINE void read_async(T* const data) const
	{
		alg_gpu_buffer<T>::read_async(data, size);
	}
	ALG_INLINE void read_sync(T* const data) const
	{
		alg_gpu_buffer<T>::read_sync(data, size);
	}
};  // struct alg_gpu_sized_buffer<T>

// ------------------------------------------------------------------------------------ //
// -----                               GPU kernels.                               ----- //
// ------------------------------------------------------------------------------------ //

template<typename U>
struct alg_opencl_kernel
{
private:
	cl::Program m_program;
	cl::Kernel m_kernel;

public:
	alg_opencl_kernel(const char* const source)
	{
		m_program = cl::Program(alg_opencl.context, {
				source, strlen(source)
		});
		if (m_program.build({ alg_opencl.default_device }) != CL_SUCCESS) {
			std::cout << " Error building: " << m_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(alg_opencl.default_device) << "\n";
			exit(1);
		}

		m_kernel = cl::Kernel(m_program, "exec");
	}

	template<typename T0>
	void exec(const size_t size, cl::Buffer& x, const T0 alpha)
	{
		ALG_GPU_CALL(m_kernel.setArg(0, x));
		ALG_GPU_CALL(m_kernel.setArg(1, alpha));
		ALG_GPU_CALL(m_kernel.setArg(2, (cl_uint)size));
		ALG_GPU_CALL(alg_opencl.queue.enqueueNDRangeKernel(m_kernel, 
			cl::NullRange, cl::NDRange(size)));
		//	cl::NullRange, cl::NDRange(alg_divide_ceil_mul(size, ALG_GPU_LOCAL_BLOCK_SIZE)), cl::NDRange(ALG_GPU_LOCAL_BLOCK_SIZE)));
	}

	template<typename T0>
	void exec(const size_t size, cl::Buffer& x, const T0 alpha, const cl::Buffer& y)
	{
		ALG_GPU_CALL(m_kernel.setArg(0, x));
		ALG_GPU_CALL(m_kernel.setArg(1, alpha));
		ALG_GPU_CALL(m_kernel.setArg(2, y));
		ALG_GPU_CALL(m_kernel.setArg(3, (cl_uint)size));
		ALG_GPU_CALL(alg_opencl.queue.enqueueNDRangeKernel(m_kernel,
			cl::NullRange, cl::NDRange(size)));
	}

	template<typename T0, typename T1>
	void exec(const size_t size, cl::Buffer& x, const T0 alpha, const cl::Buffer& y, const T1 beta, const cl::Buffer& z)
	{
		ALG_GPU_CALL(m_kernel.setArg(0, x));
		ALG_GPU_CALL(m_kernel.setArg(1, alpha));
		ALG_GPU_CALL(m_kernel.setArg(2, y));
		ALG_GPU_CALL(m_kernel.setArg(3, beta));
		ALG_GPU_CALL(m_kernel.setArg(4, z));
		ALG_GPU_CALL(m_kernel.setArg(5, (cl_uint)size));
		ALG_GPU_CALL(alg_opencl.queue.enqueueNDRangeKernel(m_kernel,
			cl::NullRange, cl::NDRange(size)));
	}
};

// ------------------------------------------------------------------------------------ //
// -----                         Vector-Vector operations.                        ----- //
// ------------------------------------------------------------------------------------ //

/// @f$ \vec x := \alpha \vec y + \beta \vec z @f$
template<typename T>
static void alg_add(alg_gpu_sized_buffer<T>& x,
	const T alpha, const alg_gpu_buffer<T>& y,
	const T beta,  const alg_gpu_buffer<T>& z)
{
	static alg_opencl_kernel<T> add_kernel(R"(
		kernel void exec(global float* const x,
			const float alpha, global const float* const y,
			const float beta,  global const float* const z,
			const unsigned size)
		{
			const size_t i = get_global_id(0);
			x[i] = alpha * y[i] + beta * z[i];
		}
	)");

	add_kernel.exec(x.size, x, alpha, y, beta, z);
}

/// @f$ \vec x += \alpha \vec y @f$
template<typename T>
static void alg_add_assign(alg_gpu_sized_buffer<T>& x,
    const T alpha, const alg_gpu_buffer<T>& y)
{
	static alg_opencl_kernel<T> add_assign_kernel(R"(
		kernel void exec(global float* const x,
			const float alpha, global const float* const y,
			const unsigned size)
		{
			const size_t i = get_global_id(0);
			//if (i >= size) return;
			x[i] += alpha * y[i];
		}
	)");

	add_assign_kernel.exec(x.size, x, alpha, y);
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/// @f$ \vec x := \alpha \vec y @f$
template<typename T>
static void alg_mul(alg_gpu_sized_buffer<T>& x,
    const T alpha, const alg_gpu_buffer<T>& y)
{
	static alg_opencl_kernel<T> mul_kernel(R"(
		kernel void exec(global float* const x,
			const float alpha, global const float* const y,
			const unsigned size)
		{
			const size_t i = get_global_id(0);
			if (i >= size) return;
			x[i] = alpha * y[i];
		}
	)");

	mul_kernel.exec(x.size, x, alpha, y);
}

/// @f$ \vec x := \alpha \vec x @f$
template<typename T>
static void alg_mul_assign(alg_gpu_sized_buffer<T>& x,
   const T alpha)
{
	static alg_opencl_kernel<T> mul_assign_kernel(R"(
		kernel void exec(global float* const x,
			const float alpha,
			const unsigned size)
		{
			const size_t i = get_global_id(0);
			if (i >= size) return;
			x[i] *= alpha;
		}
	)");

	mul_assign_kernel.exec(x.size, x, alpha);
}

// ------------------------------------------------------------------------------------ //
// -----                        Vector reduction operations.                      ----- //
// ------------------------------------------------------------------------------------ //

/// @f$ r := \sum_{i=1}^{n} y_i @f$
template<typename T>
static alg_gpu_buffer<T> alg_reduce_add(const alg_gpu_buffer<T>& y, 
	const size_t size, const size_t block_size = ALG_GPU_REDUCE_BLOCK_SIZE)
{
	static alg_opencl_kernel<T> sum_kernel(R"(
		kernel void exec(global float* const r,
			unsigned block_size, global float* const y,
			const unsigned size)
		{
			const size_t i = get_global_id(0);
			if (i == size - 1) block_size = size % 4;
			size_t j;

			r[i] = 0.0;
			for (j = i; j < block_size; ++j)
			{
				r[i] += y[j];
			}
		}
	)");

	const size_t num_blocks = alg_divide_ceil(size, block_size);
	alg_gpu_buffer<T> r(num_blocks);
	sum_kernel.exec(num_blocks, r, (cl_uint)block_size, y);

	if (num_blocks == 1)
	{
		return r;
	}
	return alg_reduce_add<T>(r, num_blocks, block_size);
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#if 0
/// @f$ r := (\vec x, \vec y) @f$
template<typename T>
static T alg_dot(const alg_gpu_sized_buffer<T>& x,
    const alg_gpu_buffer<T>& y)
{
	// https://www.youtube.com/watch?v=RzPDlnZhxtQ
	// http://developer.amd.com/resources/articles-whitepapers/opencl-optimization-case-study-simple-reductions/
	/*static alg_opencl_kernel<T> mul_kernel(R"(
		kernel void exec(
			global float* const z,
			const float _0, global float* const x,
			const float _1, global float* const y,
			const unsigned size)
		{
			const size_t i = get_global_id(0);
			if (i >= size) return;
			z[i] = x[i] * y[i];
		}
	)");

	alg_gpu_buffer<T> z(x.size);
	mul_kernel.exec(x.size, z, T(), x, T(), y);

	T r;
	const alg_gpu_buffer<T> rb = alg_reduce_add<T>(z, x.size);
	rb.read_sync(&r, 1);
	return r;*/

	static alg_opencl_kernel<T> dot_kernel(R"(
		kernel void exec(
			global float* const z,
			const float _0, global float* const x,
			const float _1, global float* const y,
			const unsigned size)
		{
			size_t i = 0;
			z[0] = 0;
			for(i = 0; i < size; ++i)
				z[0] += x[i] * y[i];
		}
	)");

	T r;
	const alg_gpu_buffer<T> rb(1);
	dot_kernel.exec(x.size, z, T(), x, T(), y);
}
#endif

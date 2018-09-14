// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
// XX                                                                                 XX //
// XX               SolverUtils - library for C++/CUDA solver development             XX //
// XX                                  Version 00.01                                  XX //
// XX                                                                                 XX //
// XX Copyright (c) 2017 Butakov Oleg                                                 XX //
// XX                                                                                 XX //
// XX Permission is hereby granted, free of charge, to any person obtaining a copy    XX //
// XX of this software and associated documentation files(the "Software"), to deal    XX //
// XX in the Software without restriction, including without limitation the rights    XX //
// XX to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       XX //
// XX copies of the Software, and to permit persons to whom the Software is           XX //
// XX furnished to do so, subject to the following conditions:                        XX //
// XX                                                                                 XX //
// XX The above copyright notice and this permission notice shall be included in all  XX //
// XX copies or substantial portions of the Software.                                 XX //
// XX                                                                                 XX //
// XX THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      XX //
// XX IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        XX //
// XX FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE      XX //
// XX AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          XX //
// XX LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   XX //
// XX OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   XX //
// XX SOFTWARE.                                                                       XX //
// XX                                                                                 XX //
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //

#pragma once

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
// XX                                  Target platform                                XX //
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //

/** 
 * We are treating non-Windows platforms as POSIX. Not exactly, but is works
 * in almost all cases. */
#if defined(_WIN32)
#	define SU_PLATFORM_WINDOWS	1
#	define SU_PLATFORM_POSIX	0
#else	// if defined(_WIN32)
#	define SU_PLATFORM_WINDOWS	0
#	define SU_PLATFORM_POSIX	1
#endif	// if defined(_WIN32)

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
// XX                                  Target compiler                                XX //
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //

/**
 * We need to recognize the Intel Compiler first - it defines MSVC's unique 
 * macros. */
#if defined(__INTEL_COMPILER)
#	define SU_COMPILER_MSC		0
#	define SU_COMPILER_GCC		0
#	define SU_COMPILER_INTEL	1
#else	// if defined(__INTEL_COMPILER)
#	define SU_COMPILER_INTEL	0
#endif	// if defined(__INTEL_COMPILER)

/**
 * As we've done before, we refer all non-MSVC compilers as GCC. */
#if defined(_MSC_VER) && !SU_COMPILER_INTEL
#	define SU_COMPILER_MSC	1
#	define SU_COMPILER_GCC	0
#else	// if defined(_MSC_VER) && !SU_COMPILER_INTEL
#	define SU_COMPILER_MSC	0
#	define SU_COMPILER_GCC	1
#endif	// if defined(_MSC_VER) && !SU_COMPILER_INTEL

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
// XX                                Target architecture                              XX //
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //

/**
 * Are we compiling solver with CUDA support? */
#if defined(__CUDACC__)
#	if !SU_COMPILER_MSC && !SU_COMPILER_GCC
#		error SovlerUtils: CUDA compiler uses MSVC or something like GCC, what is going on?
#	endif	// if !SU_COMPILER_MSC && !SU_COMPILER_GCC
#	define SU_ARCH_CUDA			1
#else	// if defined(__CUDACC__)
#	define SU_ARCH_CUDA			0
#endif	// if defined(__CUDACC__)

/**
 * If we use CUDA, are we compiling host or device code? */
#if defined(__CUDA_ARCH__)
#	define SU_ARCH_CUDA_DEVICE	1
#	define SU_ARCH_CUDA_HOST	0
#else	// if defined(__CUDACC__)
#	define SU_ARCH_CUDA_DEVICE	0
#	define SU_ARCH_CUDA_HOST	1
#endif	// if defined(__CUDACC__)

/**
 * If we are running on the CUDA device, we may have SSE support. */
#if SU_ARCH_CUDA_HOST
#	if !defined(SU_ARCH_SSE)
#		define SU_ARCH_SSE		1
#	endif	// if !defined(SU_ARCH_AVX)
#	if !defined(SU_ARCH_AVX)
#		define SU_ARCH_AVX		1
#	endif	// if !defined(SU_ARCH_AVX)
#	if !defined(SU_ARCH_AVX512)
#		define SU_ARCH_AVX512	0
#	endif	// if !defined(SU_ARCH_AVX512)
#	if SU_ARCH_AVX512 
#		if !(SU_ARCH_SSE && SU_ARCH_AVX)
#			error SovlerUtils: SSE and AVX should be enabled for AVX512 support.
#		endif	// if !(SU_ARCH_SSE && SU_ARCH_AVX)
#	endif	// if SU_ARCH_AVX512 
#	if SU_ARCH_AVX 
#		if !SU_ARCH_SSE
#			error SovlerUtils: SSE be enabled for AVX support.
#		endif	// if !SU_ARCH_SSE
#	endif	// if SU_ARCH_AVX 
#endif	// if SU_ARCH_CUDA_HOST

/**
 * Currently I don't have a machine to test AVX512 on, so let Visual Studio parse code
 * as if it is enable or supported. */
#if defined(__INTELLISENSE__) && !SU_ARCH_AVX512
#	undef SU_ARCH_AVX512
#	define SU_ARCH_AVX512	1
#endif	// if defined(__INTELLISENSE__)

/**
 * If we are running on the CUDA device, we have no SSE support. */
#if SU_ARCH_CUDA_DEVICE
#	undef SU_ARCH_SSE
#	undef SU_ARCH_AVX
#	undef SU_ARCH_AVX512
#	define SU_ARCH_SSE			0
#	define SU_ARCH_AVX			0
#	define SU_ARCH_AVX512		0
#endif	// if SU_ARCH_CUDA_DEVICE

/**
 * If we are running on the CUDA device, we may have OpenMP support. */
#if SU_ARCH_CUDA_HOST
#	if defined(_OPENMP)
#		define SU_ARCH_OPENMP	1
#	else
#		define SU_ARCH_OPENMP	0
#	endif	// if defined(_OPENMP)
#endif	// if SU_ARCH_CUDA_HOST

/**
 * If we are running on the CUDA device, we have no OpenMP support. */
#if SU_ARCH_CUDA_DEVICE
#	define SU_ARCH_OPENMP		0
#endif

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
// XX                      Code attributes & Compiler Configuration                   XX //
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //

/**
 * Configuring MSVC:
 * Enabling math constants and reducing size of Windows.h header. */
#if SU_PLATFORM_WINDOWS
#	define _USE_MATH_DEFINES	1
#	define WIN32_LEAN_AND_MEAN	1
#	define VC_EXTRALEAN			1
#endif	// if SU_PLATFORM_WINDOWS

/**
 * Wrapping compiler attributes. */
#if SU_COMPILER_MSC || SU_COMPILER_INTEL && SU_PLATFORM_WINDOWS
#	define SU_VECTORCALL		__vectorcall
#	define SU_FORCEINLINE		__forceinline
#	define SU_ALIGN_MSC(n)		__declspec(align(n))
#	define SU_ALIGN_GCC(n)
#else	// if SU_COMPILER_MSC || SU_COMPILER_INTEL && SU_PLATFORM_WINDOWS
#	define SU_VECTORCALL		
#	define SU_FORCEINLINE		__attribute__((always_inline))
#	define SU_ALIGN_MSC(n)		
#	define SU_ALIGN_GCC(n)		__attribute__((aligned(n)))
#endif	// if SU_COMPILER_MSC || SU_COMPILER_INTEL && SU_PLATFORM_WINDOWS

/**
 * Wrapping CUDA attributes. */
#if SU_ARCH_CUDA
#	define SU_CUDA_GLOBAL		__global__
#	define SU_CUDA_DEVICE		__device__
#	define SU_CUDA_HOST			__host__
#else	// if SU_ARCH_CUDA
#	define	SU_CUDA_GLOBAL
#	define	SU_CUDA_DEVICE
#	define	SU_CUDA_HOST
#endif	// if SU_ARCH_CUDA

/**
 * We often have to define both host and device functions. */
#define SU_CUDA_HOST_DEVICE SU_CUDA_HOST SU_CUDA_DEVICE

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
// XX                             Common includes and macros                          XX //
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdarg.h>

#include <float.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include <vector>
#include <algorithm>
#include <type_traits>

/**
 * CUDA includes. */
#if SU_ARCH_CUDA
#	include <cuda_runtime.h>
#	include <device_launch_parameters.h>
#endif	// if SU_ARCH_CUDA

/**
 * SSE/AVX includes. */
#if SU_ARCH_SSE
#	include <smmintrin.h>
#endif	// if SU_ARCH_SSE
#if SU_ARCH_AVX || SU_ARCH_AVX512
#	include <immintrin.h>
#endif	// if SU_ARCH_AVX || SU_ARCH_AVX512

/**
 * OpenMP includes. */
#if SU_ARCH_OPENMP
#	include <omp.h>
#endif	// if SU_ARCH_OPENMP

/** 
 * Assertions. */
#define SU_ASSERT(...) assert(__VA_ARGS__)
#define SU_ASSERT_NOT_IMPLEMENTED() SU_ASSERT(!"Not implemented")
#define SU_STATIC_ASSERT(...) static_assert(__VA_ARGS__, #__VA_ARGS__)

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //
// XX                               Library starts here...                            XX //
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX //

namespace SolverUtils
{
	namespace sustd { using namespace std; }

	using float64_t = double;
	using float32_t = float;
	using real_t = float64_t;

}	// namespace SolverUtils

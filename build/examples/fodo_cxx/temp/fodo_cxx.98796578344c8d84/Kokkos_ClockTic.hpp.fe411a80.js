var data = {lines:[
{"lineNum":"    1","line":"/*"},
{"lineNum":"    2","line":"//@HEADER"},
{"lineNum":"    3","line":"// ************************************************************************"},
{"lineNum":"    4","line":"//"},
{"lineNum":"    5","line":"//                        Kokkos v. 3.0"},
{"lineNum":"    6","line":"//       Copyright (2020) National Technology & Engineering"},
{"lineNum":"    7","line":"//               Solutions of Sandia, LLC (NTESS)."},
{"lineNum":"    8","line":"//"},
{"lineNum":"    9","line":"// Under the terms of Contract DE-NA0003525 with NTESS,"},
{"lineNum":"   10","line":"// the U.S. Government retains certain rights in this software."},
{"lineNum":"   11","line":"//"},
{"lineNum":"   12","line":"// Redistribution and use in source and binary forms, with or without"},
{"lineNum":"   13","line":"// modification, are permitted provided that the following conditions are"},
{"lineNum":"   14","line":"// met:"},
{"lineNum":"   15","line":"//"},
{"lineNum":"   16","line":"// 1. Redistributions of source code must retain the above copyright"},
{"lineNum":"   17","line":"// notice, this list of conditions and the following disclaimer."},
{"lineNum":"   18","line":"//"},
{"lineNum":"   19","line":"// 2. Redistributions in binary form must reproduce the above copyright"},
{"lineNum":"   20","line":"// notice, this list of conditions and the following disclaimer in the"},
{"lineNum":"   21","line":"// documentation and/or other materials provided with the distribution."},
{"lineNum":"   22","line":"//"},
{"lineNum":"   23","line":"// 3. Neither the name of the Corporation nor the names of the"},
{"lineNum":"   24","line":"// contributors may be used to endorse or promote products derived from"},
{"lineNum":"   25","line":"// this software without specific prior written permission."},
{"lineNum":"   26","line":"//"},
{"lineNum":"   27","line":"// THIS SOFTWARE IS PROVIDED BY NTESS \"AS IS\" AND ANY"},
{"lineNum":"   28","line":"// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE"},
{"lineNum":"   29","line":"// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR"},
{"lineNum":"   30","line":"// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE"},
{"lineNum":"   31","line":"// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,"},
{"lineNum":"   32","line":"// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,"},
{"lineNum":"   33","line":"// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR"},
{"lineNum":"   34","line":"// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF"},
{"lineNum":"   35","line":"// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING"},
{"lineNum":"   36","line":"// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS"},
{"lineNum":"   37","line":"// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."},
{"lineNum":"   38","line":"//"},
{"lineNum":"   39","line":"// Questions? Contact Christian R. Trott (crtrott@sandia.gov)"},
{"lineNum":"   40","line":"//"},
{"lineNum":"   41","line":"// ************************************************************************"},
{"lineNum":"   42","line":"//@HEADER"},
{"lineNum":"   43","line":"*/"},
{"lineNum":"   44","line":""},
{"lineNum":"   45","line":"#ifndef KOKKOS_CLOCKTIC_HPP"},
{"lineNum":"   46","line":"#define KOKKOS_CLOCKTIC_HPP"},
{"lineNum":"   47","line":""},
{"lineNum":"   48","line":"#include <Kokkos_Macros.hpp>"},
{"lineNum":"   49","line":"#include <stdint.h>"},
{"lineNum":"   50","line":"#include <chrono>"},
{"lineNum":"   51","line":"#ifdef KOKKOS_ENABLE_OPENMPTARGET"},
{"lineNum":"   52","line":"#include <omp.h>"},
{"lineNum":"   53","line":"#endif"},
{"lineNum":"   54","line":""},
{"lineNum":"   55","line":"// To use OpenCL(TM) built-in intrinsics inside kernels, we have to"},
{"lineNum":"   56","line":"// forward-declare their prototype, also see"},
{"lineNum":"   57","line":"// https://github.com/intel/pti-gpu/blob/master/chapters/binary_instrumentation/OpenCLBuiltIn.md"},
{"lineNum":"   58","line":"#if defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_ARCH_INTEL_GEN) && \\"},
{"lineNum":"   59","line":"    defined(__SYCL_DEVICE_ONLY__)"},
{"lineNum":"   60","line":"extern SYCL_EXTERNAL unsigned long __attribute__((overloadable))"},
{"lineNum":"   61","line":"intel_get_cycle_counter();"},
{"lineNum":"   62","line":"#endif"},
{"lineNum":"   63","line":""},
{"lineNum":"   64","line":"namespace Kokkos {"},
{"lineNum":"   65","line":"namespace Impl {"},
{"lineNum":"   66","line":""},
{"lineNum":"   67","line":"/**\\brief  Quick query of clock register tics"},
{"lineNum":"   68","line":" *"},
{"lineNum":"   69","line":" *  Primary use case is to, with low overhead,"},
{"lineNum":"   70","line":" *  obtain a integral value that consistently varies"},
{"lineNum":"   71","line":" *  across concurrent threads of execution within"},
{"lineNum":"   72","line":" *  a parallel algorithm."},
{"lineNum":"   73","line":" *  This value is often used to \"randomly\" seed an"},
{"lineNum":"   74","line":" *  attempt to acquire an indexed resource (e.g., bit)"},
{"lineNum":"   75","line":" *  from an array of resources (e.g., bitset) such that"},
{"lineNum":"   76","line":" *  concurrent threads will have high likelihood of"},
{"lineNum":"   77","line":" *  having different index-seed values."},
{"lineNum":"   78","line":" */"},
{"lineNum":"   79","line":""},
{"lineNum":"   80","line":"KOKKOS_FORCEINLINE_FUNCTION"},
{"lineNum":"   81","line":"uint64_t clock_tic() noexcept {"},
{"lineNum":"   82","line":"#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)"},
{"lineNum":"   83","line":""},
{"lineNum":"   84","line":"  // Return value of 64-bit hi-res clock register."},
{"lineNum":"   85","line":""},
{"lineNum":"   86","line":"  return clock64();"},
{"lineNum":"   87","line":""},
{"lineNum":"   88","line":"#elif defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_ARCH_INTEL_GEN) && \\"},
{"lineNum":"   89","line":"    defined(__SYCL_DEVICE_ONLY__)"},
{"lineNum":"   90","line":"  return intel_get_cycle_counter();"},
{"lineNum":"   91","line":"#elif defined(KOKKOS_ENABLE_OPENMPTARGET)"},
{"lineNum":"   92","line":"  return uint64_t(omp_get_wtime() * 1.e9);"},
{"lineNum":"   93","line":"#elif defined(__i386__) || defined(__x86_64)"},
{"lineNum":"   94","line":""},
{"lineNum":"   95","line":"  // Return value of 64-bit hi-res clock register."},
{"lineNum":"   96","line":""},
{"lineNum":"   97","line":"  unsigned a = 0, d = 0;"},
{"lineNum":"   98","line":""},
{"lineNum":"   99","line":"  __asm__ volatile(\"rdtsc\" : \"=a\"(a), \"=d\"(d));","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  100","line":""},
{"lineNum":"  101","line":"  return ((uint64_t)a) | (((uint64_t)d) << 32);"},
{"lineNum":"  102","line":""},
{"lineNum":"  103","line":"#elif defined(__powerpc) || defined(__powerpc__) || defined(__powerpc64__) || \\"},
{"lineNum":"  104","line":"    defined(__POWERPC__) || defined(__ppc__) || defined(__ppc64__)"},
{"lineNum":"  105","line":""},
{"lineNum":"  106","line":"  unsigned int cycles = 0;"},
{"lineNum":"  107","line":""},
{"lineNum":"  108","line":"  asm volatile(\"mftb %0\" : \"=r\"(cycles));"},
{"lineNum":"  109","line":""},
{"lineNum":"  110","line":"  return (uint64_t)cycles;"},
{"lineNum":"  111","line":""},
{"lineNum":"  112","line":"#else"},
{"lineNum":"  113","line":""},
{"lineNum":"  114","line":"  return (uint64_t)std::chrono::high_resolution_clock::now()"},
{"lineNum":"  115","line":"      .time_since_epoch()"},
{"lineNum":"  116","line":"      .count();"},
{"lineNum":"  117","line":""},
{"lineNum":"  118","line":"#endif"},
{"lineNum":"  119","line":"}"},
{"lineNum":"  120","line":""},
{"lineNum":"  121","line":"}  // namespace Impl"},
{"lineNum":"  122","line":"}  // namespace Kokkos"},
{"lineNum":"  123","line":""},
{"lineNum":"  124","line":"#endif  // KOKKOS_CLOCKTIC_HPP"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 1, "covered" : 0,};
var merged_data = [];

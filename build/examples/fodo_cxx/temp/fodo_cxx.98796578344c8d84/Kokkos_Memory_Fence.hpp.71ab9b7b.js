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
{"lineNum":"   39","line":"// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)"},
{"lineNum":"   40","line":"//"},
{"lineNum":"   41","line":"// ************************************************************************"},
{"lineNum":"   42","line":"//@HEADER"},
{"lineNum":"   43","line":"*/"},
{"lineNum":"   44","line":""},
{"lineNum":"   45","line":"#include <Kokkos_Macros.hpp>"},
{"lineNum":"   46","line":"#if defined(KOKKOS_ATOMIC_HPP) && !defined(KOKKOS_MEMORY_FENCE_HPP)"},
{"lineNum":"   47","line":"#define KOKKOS_MEMORY_FENCE_HPP"},
{"lineNum":"   48","line":"namespace Kokkos {"},
{"lineNum":"   49","line":""},
{"lineNum":"   50","line":"//----------------------------------------------------------------------------"},
{"lineNum":"   51","line":""},
{"lineNum":"   52","line":"KOKKOS_FORCEINLINE_FUNCTION"},
{"lineNum":"   53","line":"void memory_fence() {"},
{"lineNum":"   54","line":"#if defined(__CUDA_ARCH__)"},
{"lineNum":"   55","line":"  __threadfence();"},
{"lineNum":"   56","line":"#elif defined(KOKKOS_ENABLE_OPENMPTARGET)"},
{"lineNum":"   57","line":"#pragma omp flush"},
{"lineNum":"   58","line":"#elif defined(__HIP_DEVICE_COMPILE__)"},
{"lineNum":"   59","line":"  __threadfence();"},
{"lineNum":"   60","line":"#elif defined(KOKKOS_ENABLE_SYCL) && defined(__SYCL_DEVICE_ONLY__)"},
{"lineNum":"   61","line":"  sycl::ONEAPI::atomic_fence(sycl::ONEAPI::memory_order::acq_rel,"},
{"lineNum":"   62","line":"                             sycl::ONEAPI::memory_scope::device);"},
{"lineNum":"   63","line":"#elif defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)"},
{"lineNum":"   64","line":"  asm volatile(\"mfence\" ::: \"memory\");"},
{"lineNum":"   65","line":"#elif defined(KOKKOS_ENABLE_GNU_ATOMICS) || \\"},
{"lineNum":"   66","line":"    (defined(KOKKOS_COMPILER_NVCC) && defined(KOKKOS_ENABLE_INTEL_ATOMICS))"},
{"lineNum":"   67","line":"  __sync_synchronize();","class":"lineNoCov","hits":"0","possible_hits":"36",},
{"lineNum":"   68","line":"#elif defined(KOKKOS_ENABLE_INTEL_ATOMICS)"},
{"lineNum":"   69","line":"  _mm_mfence();"},
{"lineNum":"   70","line":"#elif defined(KOKKOS_ENABLE_OPENMP_ATOMICS)"},
{"lineNum":"   71","line":"#pragma omp flush"},
{"lineNum":"   72","line":"#elif defined(KOKKOS_ENABLE_WINDOWS_ATOMICS)"},
{"lineNum":"   73","line":"  MemoryBarrier();"},
{"lineNum":"   74","line":"#elif !defined(KOKKOS_ENABLE_SERIAL_ATOMICS)"},
{"lineNum":"   75","line":"#error \"Error: memory_fence() not defined\""},
{"lineNum":"   76","line":"#endif"},
{"lineNum":"   77","line":"}"},
{"lineNum":"   78","line":""},
{"lineNum":"   79","line":"//////////////////////////////////////////////////////"},
{"lineNum":"   80","line":"// store_fence()"},
{"lineNum":"   81","line":"//"},
{"lineNum":"   82","line":"// If possible use a store fence on the architecture, if not run a full memory"},
{"lineNum":"   83","line":"// fence"},
{"lineNum":"   84","line":""},
{"lineNum":"   85","line":"KOKKOS_FORCEINLINE_FUNCTION"},
{"lineNum":"   86","line":"void store_fence() {"},
{"lineNum":"   87","line":"#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)"},
{"lineNum":"   88","line":"  asm volatile(\"sfence\" ::: \"memory\");"},
{"lineNum":"   89","line":"#else"},
{"lineNum":"   90","line":"  memory_fence();"},
{"lineNum":"   91","line":"#endif"},
{"lineNum":"   92","line":"}"},
{"lineNum":"   93","line":""},
{"lineNum":"   94","line":"//////////////////////////////////////////////////////"},
{"lineNum":"   95","line":"// load_fence()"},
{"lineNum":"   96","line":"//"},
{"lineNum":"   97","line":"// If possible use a load fence on the architecture, if not run a full memory"},
{"lineNum":"   98","line":"// fence"},
{"lineNum":"   99","line":""},
{"lineNum":"  100","line":"KOKKOS_FORCEINLINE_FUNCTION"},
{"lineNum":"  101","line":"void load_fence() {"},
{"lineNum":"  102","line":"#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)"},
{"lineNum":"  103","line":"  asm volatile(\"lfence\" ::: \"memory\");"},
{"lineNum":"  104","line":"#else"},
{"lineNum":"  105","line":"  memory_fence();"},
{"lineNum":"  106","line":"#endif"},
{"lineNum":"  107","line":"}"},
{"lineNum":"  108","line":""},
{"lineNum":"  109","line":"}  // namespace Kokkos"},
{"lineNum":"  110","line":""},
{"lineNum":"  111","line":"#endif"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 1, "covered" : 0,};
var merged_data = [];

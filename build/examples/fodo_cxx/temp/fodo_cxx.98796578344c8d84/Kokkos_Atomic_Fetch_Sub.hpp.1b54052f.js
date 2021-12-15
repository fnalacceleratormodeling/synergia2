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
{"lineNum":"   45","line":"#if defined(KOKKOS_ENABLE_RFO_PREFETCH)"},
{"lineNum":"   46","line":"#include <xmmintrin.h>"},
{"lineNum":"   47","line":"#endif"},
{"lineNum":"   48","line":""},
{"lineNum":"   49","line":"#include <Kokkos_Macros.hpp>"},
{"lineNum":"   50","line":"#if defined(KOKKOS_ATOMIC_HPP) && !defined(KOKKOS_ATOMIC_FETCH_SUB_HPP)"},
{"lineNum":"   51","line":"#define KOKKOS_ATOMIC_FETCH_SUB_HPP"},
{"lineNum":"   52","line":""},
{"lineNum":"   53","line":"#if defined(KOKKOS_ENABLE_CUDA)"},
{"lineNum":"   54","line":"#include <Cuda/Kokkos_Cuda_Version_9_8_Compatibility.hpp>"},
{"lineNum":"   55","line":"#endif"},
{"lineNum":"   56","line":""},
{"lineNum":"   57","line":"namespace Kokkos {"},
{"lineNum":"   58","line":""},
{"lineNum":"   59","line":"//----------------------------------------------------------------------------"},
{"lineNum":"   60","line":""},
{"lineNum":"   61","line":"#if defined(KOKKOS_ENABLE_CUDA)"},
{"lineNum":"   62","line":"#if defined(__CUDA_ARCH__) || defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND)"},
{"lineNum":"   63","line":""},
{"lineNum":"   64","line":"// Support for int, unsigned int, unsigned long long int, and float"},
{"lineNum":"   65","line":""},
{"lineNum":"   66","line":"__inline__ __device__ int atomic_fetch_sub(volatile int* const dest,"},
{"lineNum":"   67","line":"                                           const int val) {"},
{"lineNum":"   68","line":"  return atomicSub((int*)dest, val);"},
{"lineNum":"   69","line":"}"},
{"lineNum":"   70","line":""},
{"lineNum":"   71","line":"__inline__ __device__ unsigned int atomic_fetch_sub("},
{"lineNum":"   72","line":"    volatile unsigned int* const dest, const unsigned int val) {"},
{"lineNum":"   73","line":"  return atomicSub((unsigned int*)dest, val);"},
{"lineNum":"   74","line":"}"},
{"lineNum":"   75","line":""},
{"lineNum":"   76","line":"__inline__ __device__ unsigned int atomic_fetch_sub("},
{"lineNum":"   77","line":"    volatile int64_t* const dest, const int64_t val) {"},
{"lineNum":"   78","line":"  return atomic_fetch_add(dest, -val);"},
{"lineNum":"   79","line":"}"},
{"lineNum":"   80","line":""},
{"lineNum":"   81","line":"__inline__ __device__ unsigned int atomic_fetch_sub(volatile float* const dest,"},
{"lineNum":"   82","line":"                                                    const float val) {"},
{"lineNum":"   83","line":"  return atomicAdd((float*)dest, -val);"},
{"lineNum":"   84","line":"}"},
{"lineNum":"   85","line":""},
{"lineNum":"   86","line":"#if (600 <= __CUDA_ARCH__)"},
{"lineNum":"   87","line":"__inline__ __device__ unsigned int atomic_fetch_sub(volatile double* const dest,"},
{"lineNum":"   88","line":"                                                    const double val) {"},
{"lineNum":"   89","line":"  return atomicAdd((double*)dest, -val);"},
{"lineNum":"   90","line":"}"},
{"lineNum":"   91","line":"#endif"},
{"lineNum":"   92","line":""},
{"lineNum":"   93","line":"template <typename T>"},
{"lineNum":"   94","line":"__inline__ __device__ T atomic_fetch_sub("},
{"lineNum":"   95","line":"    volatile T* const dest,"},
{"lineNum":"   96","line":"    typename std::enable_if<sizeof(T) == sizeof(int), const T>::type val) {"},
{"lineNum":"   97","line":"  union U {"},
{"lineNum":"   98","line":"    int i;"},
{"lineNum":"   99","line":"    T t;"},
{"lineNum":"  100","line":"    KOKKOS_INLINE_FUNCTION U() {}"},
{"lineNum":"  101","line":"  } oldval, assume, newval;"},
{"lineNum":"  102","line":""},
{"lineNum":"  103","line":"  oldval.t = *dest;"},
{"lineNum":"  104","line":""},
{"lineNum":"  105","line":"  do {"},
{"lineNum":"  106","line":"    assume.i = oldval.i;"},
{"lineNum":"  107","line":"    newval.t = assume.t - val;"},
{"lineNum":"  108","line":"    oldval.i = atomicCAS((int*)dest, assume.i, newval.i);"},
{"lineNum":"  109","line":"  } while (assume.i != oldval.i);"},
{"lineNum":"  110","line":""},
{"lineNum":"  111","line":"  return oldval.t;"},
{"lineNum":"  112","line":"}"},
{"lineNum":"  113","line":""},
{"lineNum":"  114","line":"template <typename T>"},
{"lineNum":"  115","line":"__inline__ __device__ T atomic_fetch_sub("},
{"lineNum":"  116","line":"    volatile T* const dest,"},
{"lineNum":"  117","line":"    typename std::enable_if<sizeof(T) != sizeof(int) &&"},
{"lineNum":"  118","line":"                                sizeof(T) == sizeof(unsigned long long int),"},
{"lineNum":"  119","line":"                            const T>::type val) {"},
{"lineNum":"  120","line":"  union U {"},
{"lineNum":"  121","line":"    unsigned long long int i;"},
{"lineNum":"  122","line":"    T t;"},
{"lineNum":"  123","line":"    KOKKOS_INLINE_FUNCTION U() {}"},
{"lineNum":"  124","line":"  } oldval, assume, newval;"},
{"lineNum":"  125","line":""},
{"lineNum":"  126","line":"  oldval.t = *dest;"},
{"lineNum":"  127","line":""},
{"lineNum":"  128","line":"  do {"},
{"lineNum":"  129","line":"    assume.i = oldval.i;"},
{"lineNum":"  130","line":"    newval.t = assume.t - val;"},
{"lineNum":"  131","line":"    oldval.i = atomicCAS((unsigned long long int*)dest, assume.i, newval.i);"},
{"lineNum":"  132","line":"  } while (assume.i != oldval.i);"},
{"lineNum":"  133","line":""},
{"lineNum":"  134","line":"  return oldval.t;"},
{"lineNum":"  135","line":"}"},
{"lineNum":"  136","line":""},
{"lineNum":"  137","line":"//----------------------------------------------------------------------------"},
{"lineNum":"  138","line":""},
{"lineNum":"  139","line":"template <typename T>"},
{"lineNum":"  140","line":"__inline__ __device__ T"},
{"lineNum":"  141","line":"atomic_fetch_sub(volatile T* const dest,"},
{"lineNum":"  142","line":"                 typename std::enable_if<(sizeof(T) != 4) && (sizeof(T) != 8),"},
{"lineNum":"  143","line":"                                         const T>::type& val) {"},
{"lineNum":"  144","line":"  T return_val;"},
{"lineNum":"  145","line":"  // This is a way to (hopefully) avoid dead lock in a warp"},
{"lineNum":"  146","line":"  int done = 0;"},
{"lineNum":"  147","line":"#ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK"},
{"lineNum":"  148","line":"  unsigned int mask   = KOKKOS_IMPL_CUDA_ACTIVEMASK;"},
{"lineNum":"  149","line":"  unsigned int active = KOKKOS_IMPL_CUDA_BALLOT_MASK(mask, 1);"},
{"lineNum":"  150","line":"#else"},
{"lineNum":"  151","line":"  unsigned int active = KOKKOS_IMPL_CUDA_BALLOT(1);"},
{"lineNum":"  152","line":"#endif"},
{"lineNum":"  153","line":"  unsigned int done_active = 0;"},
{"lineNum":"  154","line":"  while (active != done_active) {"},
{"lineNum":"  155","line":"    if (!done) {"},
{"lineNum":"  156","line":"      if (Impl::lock_address_cuda_space((void*)dest)) {"},
{"lineNum":"  157","line":"        Kokkos::memory_fence();"},
{"lineNum":"  158","line":"        return_val = *dest;"},
{"lineNum":"  159","line":"        *dest      = return_val - val;"},
{"lineNum":"  160","line":"        Kokkos::memory_fence();"},
{"lineNum":"  161","line":"        Impl::unlock_address_cuda_space((void*)dest);"},
{"lineNum":"  162","line":"        done = 1;"},
{"lineNum":"  163","line":"      }"},
{"lineNum":"  164","line":"    }"},
{"lineNum":"  165","line":"#ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK"},
{"lineNum":"  166","line":"    done_active = KOKKOS_IMPL_CUDA_BALLOT_MASK(mask, done);"},
{"lineNum":"  167","line":"#else"},
{"lineNum":"  168","line":"    done_active = KOKKOS_IMPL_CUDA_BALLOT(done);"},
{"lineNum":"  169","line":"#endif"},
{"lineNum":"  170","line":"  }"},
{"lineNum":"  171","line":"  return return_val;"},
{"lineNum":"  172","line":"}"},
{"lineNum":"  173","line":"#endif"},
{"lineNum":"  174","line":"#endif"},
{"lineNum":"  175","line":"//----------------------------------------------------------------------------"},
{"lineNum":"  176","line":"#if !defined(__CUDA_ARCH__) || defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND)"},
{"lineNum":"  177","line":"#if defined(KOKKOS_ENABLE_GNU_ATOMICS) || defined(KOKKOS_ENABLE_INTEL_ATOMICS)"},
{"lineNum":"  178","line":""},
{"lineNum":"  179","line":"inline int atomic_fetch_sub(volatile int* const dest, const int val) {"},
{"lineNum":"  180","line":"#if defined(KOKKOS_ENABLE_RFO_PREFETCH)"},
{"lineNum":"  181","line":"  _mm_prefetch((const char*)dest, _MM_HINT_ET0);"},
{"lineNum":"  182","line":"#endif"},
{"lineNum":"  183","line":"  return __sync_fetch_and_sub(dest, val);","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  184","line":"}"},
{"lineNum":"  185","line":""},
{"lineNum":"  186","line":"inline long int atomic_fetch_sub(volatile long int* const dest,"},
{"lineNum":"  187","line":"                                 const long int val) {"},
{"lineNum":"  188","line":"#if defined(KOKKOS_ENABLE_RFO_PREFETCH)"},
{"lineNum":"  189","line":"  _mm_prefetch((const char*)dest, _MM_HINT_ET0);"},
{"lineNum":"  190","line":"#endif"},
{"lineNum":"  191","line":"  return __sync_fetch_and_sub(dest, val);"},
{"lineNum":"  192","line":"}"},
{"lineNum":"  193","line":""},
{"lineNum":"  194","line":"#if defined(KOKKOS_ENABLE_GNU_ATOMICS)"},
{"lineNum":"  195","line":""},
{"lineNum":"  196","line":"inline unsigned int atomic_fetch_sub(volatile unsigned int* const dest,"},
{"lineNum":"  197","line":"                                     const unsigned int val) {"},
{"lineNum":"  198","line":"#if defined(KOKKOS_ENABLE_RFO_PREFETCH)"},
{"lineNum":"  199","line":"  _mm_prefetch((const char*)dest, _MM_HINT_ET0);"},
{"lineNum":"  200","line":"#endif"},
{"lineNum":"  201","line":"  return __sync_fetch_and_sub(dest, val);"},
{"lineNum":"  202","line":"}"},
{"lineNum":"  203","line":""},
{"lineNum":"  204","line":"inline unsigned long int atomic_fetch_sub("},
{"lineNum":"  205","line":"    volatile unsigned long int* const dest, const unsigned long int val) {"},
{"lineNum":"  206","line":"#if defined(KOKKOS_ENABLE_RFO_PREFETCH)"},
{"lineNum":"  207","line":"  _mm_prefetch((const char*)dest, _MM_HINT_ET0);"},
{"lineNum":"  208","line":"#endif"},
{"lineNum":"  209","line":"  return __sync_fetch_and_sub(dest, val);"},
{"lineNum":"  210","line":"}"},
{"lineNum":"  211","line":""},
{"lineNum":"  212","line":"inline unsigned long long int atomic_fetch_sub("},
{"lineNum":"  213","line":"    volatile unsigned long long int* const dest,"},
{"lineNum":"  214","line":"    const unsigned long long int val) {"},
{"lineNum":"  215","line":"#if defined(KOKKOS_ENABLE_RFO_PREFETCH)"},
{"lineNum":"  216","line":"  _mm_prefetch((const char*)dest, _MM_HINT_ET0);"},
{"lineNum":"  217","line":"#endif"},
{"lineNum":"  218","line":"  return __sync_fetch_and_sub(dest, val);"},
{"lineNum":"  219","line":"}"},
{"lineNum":"  220","line":""},
{"lineNum":"  221","line":"#endif"},
{"lineNum":"  222","line":""},
{"lineNum":"  223","line":"template <typename T>"},
{"lineNum":"  224","line":"inline T atomic_fetch_sub("},
{"lineNum":"  225","line":"    volatile T* const dest,"},
{"lineNum":"  226","line":"    typename std::enable_if<sizeof(T) == sizeof(int), const T>::type val) {"},
{"lineNum":"  227","line":"  union U {"},
{"lineNum":"  228","line":"    int i;"},
{"lineNum":"  229","line":"    T t;"},
{"lineNum":"  230","line":"    KOKKOS_INLINE_FUNCTION U() {}"},
{"lineNum":"  231","line":"  } oldval, assume, newval;"},
{"lineNum":"  232","line":""},
{"lineNum":"  233","line":"#if defined(KOKKOS_ENABLE_RFO_PREFETCH)"},
{"lineNum":"  234","line":"  _mm_prefetch((const char*)dest, _MM_HINT_ET0);"},
{"lineNum":"  235","line":"#endif"},
{"lineNum":"  236","line":""},
{"lineNum":"  237","line":"  oldval.t = *dest;"},
{"lineNum":"  238","line":""},
{"lineNum":"  239","line":"  do {"},
{"lineNum":"  240","line":"    assume.i = oldval.i;"},
{"lineNum":"  241","line":"    newval.t = assume.t - val;"},
{"lineNum":"  242","line":"    oldval.i = __sync_val_compare_and_swap((int*)dest, assume.i, newval.i);"},
{"lineNum":"  243","line":"  } while (assume.i != oldval.i);"},
{"lineNum":"  244","line":""},
{"lineNum":"  245","line":"  return oldval.t;"},
{"lineNum":"  246","line":"}"},
{"lineNum":"  247","line":""},
{"lineNum":"  248","line":"template <typename T>"},
{"lineNum":"  249","line":"inline T atomic_fetch_sub(volatile T* const dest,"},
{"lineNum":"  250","line":"                          typename std::enable_if<sizeof(T) != sizeof(int) &&"},
{"lineNum":"  251","line":"                                                      sizeof(T) == sizeof(long),"},
{"lineNum":"  252","line":"                                                  const T>::type val) {"},
{"lineNum":"  253","line":"#if defined(KOKKOS_ENABLE_RFO_PREFETCH)"},
{"lineNum":"  254","line":"  _mm_prefetch((const char*)dest, _MM_HINT_ET0);"},
{"lineNum":"  255","line":"#endif"},
{"lineNum":"  256","line":""},
{"lineNum":"  257","line":"  union U {"},
{"lineNum":"  258","line":"    long i;"},
{"lineNum":"  259","line":"    T t;"},
{"lineNum":"  260","line":"    KOKKOS_INLINE_FUNCTION U() {}"},
{"lineNum":"  261","line":"  } oldval, assume, newval;"},
{"lineNum":"  262","line":""},
{"lineNum":"  263","line":"  oldval.t = *dest;"},
{"lineNum":"  264","line":""},
{"lineNum":"  265","line":"  do {"},
{"lineNum":"  266","line":"    assume.i = oldval.i;"},
{"lineNum":"  267","line":"    newval.t = assume.t - val;"},
{"lineNum":"  268","line":"    oldval.i = __sync_val_compare_and_swap((long*)dest, assume.i, newval.i);"},
{"lineNum":"  269","line":"  } while (assume.i != oldval.i);"},
{"lineNum":"  270","line":""},
{"lineNum":"  271","line":"  return oldval.t;"},
{"lineNum":"  272","line":"}"},
{"lineNum":"  273","line":""},
{"lineNum":"  274","line":"//----------------------------------------------------------------------------"},
{"lineNum":"  275","line":""},
{"lineNum":"  276","line":"template <typename T>"},
{"lineNum":"  277","line":"inline T atomic_fetch_sub("},
{"lineNum":"  278","line":"    volatile T* const dest,"},
{"lineNum":"  279","line":"    typename std::enable_if<(sizeof(T) != 4) && (sizeof(T) != 8),"},
{"lineNum":"  280","line":"                            const T>::type& val) {"},
{"lineNum":"  281","line":"#if defined(KOKKOS_ENABLE_RFO_PREFETCH)"},
{"lineNum":"  282","line":"  _mm_prefetch((const char*)dest, _MM_HINT_ET0);"},
{"lineNum":"  283","line":"#endif"},
{"lineNum":"  284","line":""},
{"lineNum":"  285","line":"  while (!Impl::lock_address_host_space((void*)dest))"},
{"lineNum":"  286","line":"    ;"},
{"lineNum":"  287","line":"  Kokkos::memory_fence();"},
{"lineNum":"  288","line":"  T return_val = *dest;"},
{"lineNum":"  289","line":"  *dest        = return_val - val;"},
{"lineNum":"  290","line":"  Kokkos::memory_fence();"},
{"lineNum":"  291","line":"  Impl::unlock_address_host_space((void*)dest);"},
{"lineNum":"  292","line":"  return return_val;"},
{"lineNum":"  293","line":"}"},
{"lineNum":"  294","line":""},
{"lineNum":"  295","line":"//----------------------------------------------------------------------------"},
{"lineNum":"  296","line":""},
{"lineNum":"  297","line":"#elif defined(KOKKOS_ENABLE_OPENMP_ATOMICS)"},
{"lineNum":"  298","line":""},
{"lineNum":"  299","line":"template <typename T>"},
{"lineNum":"  300","line":"T atomic_fetch_sub(volatile T* const dest, const T val) {"},
{"lineNum":"  301","line":"  T retval;"},
{"lineNum":"  302","line":"#pragma omp atomic capture"},
{"lineNum":"  303","line":"  {"},
{"lineNum":"  304","line":"    retval = dest[0];"},
{"lineNum":"  305","line":"    dest[0] -= val;"},
{"lineNum":"  306","line":"  }"},
{"lineNum":"  307","line":"  return retval;"},
{"lineNum":"  308","line":"}"},
{"lineNum":"  309","line":""},
{"lineNum":"  310","line":"#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)"},
{"lineNum":"  311","line":""},
{"lineNum":"  312","line":"template <typename T>"},
{"lineNum":"  313","line":"T atomic_fetch_sub(volatile T* const dest_v, const T val) {"},
{"lineNum":"  314","line":"  T* dest  = const_cast<T*>(dest_v);"},
{"lineNum":"  315","line":"  T retval = *dest;"},
{"lineNum":"  316","line":"  *dest -= val;"},
{"lineNum":"  317","line":"  return retval;"},
{"lineNum":"  318","line":"}"},
{"lineNum":"  319","line":""},
{"lineNum":"  320","line":"#endif"},
{"lineNum":"  321","line":"#endif"},
{"lineNum":"  322","line":""},
{"lineNum":"  323","line":"// dummy for non-CUDA Kokkos headers being processed by NVCC"},
{"lineNum":"  324","line":"#if defined(__CUDA_ARCH__) && !defined(KOKKOS_ENABLE_CUDA)"},
{"lineNum":"  325","line":"template <typename T>"},
{"lineNum":"  326","line":"__inline__ __device__ T atomic_fetch_sub(volatile T* const,"},
{"lineNum":"  327","line":"                                         Kokkos::Impl::identity_t<T>) {"},
{"lineNum":"  328","line":"  return T();"},
{"lineNum":"  329","line":"}"},
{"lineNum":"  330","line":"#endif"},
{"lineNum":"  331","line":""},
{"lineNum":"  332","line":"}  // namespace Kokkos"},
{"lineNum":"  333","line":""},
{"lineNum":"  334","line":"#include <impl/Kokkos_Atomic_Assembly.hpp>"},
{"lineNum":"  335","line":"#endif"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 1, "covered" : 0,};
var merged_data = [];

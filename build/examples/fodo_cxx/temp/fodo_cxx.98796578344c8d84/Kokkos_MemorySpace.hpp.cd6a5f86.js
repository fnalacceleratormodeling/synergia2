var data = {lines:[
{"lineNum":"    1","line":"/*"},
{"lineNum":"    2","line":"//@HEADER"},
{"lineNum":"    3","line":"// ************************************************************************"},
{"lineNum":"    4","line":"//"},
{"lineNum":"    5","line":"//                        Kokkos v. 3.0"},
{"lineNum":"    6","line":"//              Copyright (2019) Sandia Corporation"},
{"lineNum":"    7","line":"//"},
{"lineNum":"    8","line":"// Under the terms of Contract DE-NA0003525 with NTESS,"},
{"lineNum":"    9","line":"// the U.S. Government retains certain rights in this software."},
{"lineNum":"   10","line":"//"},
{"lineNum":"   11","line":"// Redistribution and use in source and binary forms, with or without"},
{"lineNum":"   12","line":"// modification, are permitted provided that the following conditions are"},
{"lineNum":"   13","line":"// met:"},
{"lineNum":"   14","line":"//"},
{"lineNum":"   15","line":"// 1. Redistributions of source code must retain the above copyright"},
{"lineNum":"   16","line":"// notice, this list of conditions and the following disclaimer."},
{"lineNum":"   17","line":"//"},
{"lineNum":"   18","line":"// 2. Redistributions in binary form must reproduce the above copyright"},
{"lineNum":"   19","line":"// notice, this list of conditions and the following disclaimer in the"},
{"lineNum":"   20","line":"// documentation and/or other materials provided with the distribution."},
{"lineNum":"   21","line":"//"},
{"lineNum":"   22","line":"// 3. Neither the name of the Corporation nor the names of the"},
{"lineNum":"   23","line":"// contributors may be used to endorse or promote products derived from"},
{"lineNum":"   24","line":"// this software without specific prior written permission."},
{"lineNum":"   25","line":"//"},
{"lineNum":"   26","line":"// THIS SOFTWARE IS PROVIDED BY NTESS \"AS IS\" AND ANY"},
{"lineNum":"   27","line":"// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE"},
{"lineNum":"   28","line":"// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR"},
{"lineNum":"   29","line":"// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE"},
{"lineNum":"   30","line":"// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,"},
{"lineNum":"   31","line":"// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,"},
{"lineNum":"   32","line":"// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR"},
{"lineNum":"   33","line":"// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF"},
{"lineNum":"   34","line":"// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING"},
{"lineNum":"   35","line":"// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS"},
{"lineNum":"   36","line":"// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."},
{"lineNum":"   37","line":"//"},
{"lineNum":"   38","line":"// Questions? Contact Christian R. Trott (crtrott@sandia.gov)"},
{"lineNum":"   39","line":"//"},
{"lineNum":"   40","line":"// ************************************************************************"},
{"lineNum":"   41","line":"//@HEADER"},
{"lineNum":"   42","line":"*/"},
{"lineNum":"   43","line":""},
{"lineNum":"   44","line":"/** @file Kokkos_MemorySpace.hpp"},
{"lineNum":"   45","line":" *"},
{"lineNum":"   46","line":" *  Operations common to memory space instances, or at least default"},
{"lineNum":"   47","line":" *  implementations thereof."},
{"lineNum":"   48","line":" */"},
{"lineNum":"   49","line":""},
{"lineNum":"   50","line":"#ifndef KOKKOS_IMPL_MEMORYSPACE_HPP"},
{"lineNum":"   51","line":"#define KOKKOS_IMPL_MEMORYSPACE_HPP"},
{"lineNum":"   52","line":""},
{"lineNum":"   53","line":"#include <Kokkos_Macros.hpp>"},
{"lineNum":"   54","line":"#include <impl/Kokkos_SharedAlloc.hpp>"},
{"lineNum":"   55","line":"#include <impl/Kokkos_Error.hpp>"},
{"lineNum":"   56","line":""},
{"lineNum":"   57","line":"#include <string>"},
{"lineNum":"   58","line":""},
{"lineNum":"   59","line":"namespace Kokkos {"},
{"lineNum":"   60","line":"namespace Impl {"},
{"lineNum":"   61","line":""},
{"lineNum":"   62","line":"// Defined in implementation file to avoid having to include iostream"},
{"lineNum":"   63","line":"void safe_throw_allocation_with_header_failure("},
{"lineNum":"   64","line":"    std::string const &space_name, std::string const &label,"},
{"lineNum":"   65","line":"    Kokkos::Experimental::RawMemoryAllocationFailure const &failure);"},
{"lineNum":"   66","line":""},
{"lineNum":"   67","line":"template <class MemorySpace>"},
{"lineNum":"   68","line":"SharedAllocationHeader *checked_allocation_with_header(MemorySpace const &space,"},
{"lineNum":"   69","line":"                                                       std::string const &label,"},
{"lineNum":"   70","line":"                                                       size_t alloc_size) {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   71","line":"  try {"},
{"lineNum":"   72","line":"    return reinterpret_cast<SharedAllocationHeader *>(space.allocate("},
{"lineNum":"   73","line":"        label.c_str(), alloc_size + sizeof(SharedAllocationHeader),","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   74","line":"        alloc_size));"},
{"lineNum":"   75","line":"  } catch (Kokkos::Experimental::RawMemoryAllocationFailure const &failure) {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   76","line":"    safe_throw_allocation_with_header_failure(space.name(), label, failure);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   77","line":"  }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   78","line":"  return nullptr;  // unreachable"},
{"lineNum":"   79","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   80","line":""},
{"lineNum":"   81","line":"}  // end namespace Impl"},
{"lineNum":"   82","line":"}  // end namespace Kokkos"},
{"lineNum":"   83","line":""},
{"lineNum":"   84","line":"#endif  // KOKKOS_IMPL_MEMORYSPACE_HPP"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 6, "covered" : 0,};
var merged_data = [];

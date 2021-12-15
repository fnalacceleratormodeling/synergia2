var data = {lines:[
{"lineNum":"    1","line":"// -*- C++ -*-"},
{"lineNum":"    2","line":"//===----------------------------------------------------------------------===//"},
{"lineNum":"    3","line":"//"},
{"lineNum":"    4","line":"// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions."},
{"lineNum":"    5","line":"// See https://llvm.org/LICENSE.txt for license information."},
{"lineNum":"    6","line":"// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception"},
{"lineNum":"    7","line":"//"},
{"lineNum":"    8","line":"//===----------------------------------------------------------------------===//"},
{"lineNum":"    9","line":""},
{"lineNum":"   10","line":"#ifndef _LIBCPP___ITERATOR_DISTANCE_H"},
{"lineNum":"   11","line":"#define _LIBCPP___ITERATOR_DISTANCE_H"},
{"lineNum":"   12","line":""},
{"lineNum":"   13","line":"#include <__config>"},
{"lineNum":"   14","line":"#include <__iterator/iterator_traits.h>"},
{"lineNum":"   15","line":""},
{"lineNum":"   16","line":"#if !defined(_LIBCPP_HAS_NO_PRAGMA_SYSTEM_HEADER)"},
{"lineNum":"   17","line":"#pragma GCC system_header"},
{"lineNum":"   18","line":"#endif"},
{"lineNum":"   19","line":""},
{"lineNum":"   20","line":"_LIBCPP_PUSH_MACROS"},
{"lineNum":"   21","line":"#include <__undef_macros>"},
{"lineNum":"   22","line":""},
{"lineNum":"   23","line":"_LIBCPP_BEGIN_NAMESPACE_STD"},
{"lineNum":"   24","line":""},
{"lineNum":"   25","line":"template <class _InputIter>"},
{"lineNum":"   26","line":"inline _LIBCPP_INLINE_VISIBILITY _LIBCPP_CONSTEXPR_AFTER_CXX14"},
{"lineNum":"   27","line":"typename iterator_traits<_InputIter>::difference_type"},
{"lineNum":"   28","line":"__distance(_InputIter __first, _InputIter __last, input_iterator_tag)"},
{"lineNum":"   29","line":"{"},
{"lineNum":"   30","line":"    typename iterator_traits<_InputIter>::difference_type __r(0);"},
{"lineNum":"   31","line":"    for (; __first != __last; ++__first)"},
{"lineNum":"   32","line":"        ++__r;"},
{"lineNum":"   33","line":"    return __r;"},
{"lineNum":"   34","line":"}"},
{"lineNum":"   35","line":""},
{"lineNum":"   36","line":"template <class _RandIter>"},
{"lineNum":"   37","line":"inline _LIBCPP_INLINE_VISIBILITY _LIBCPP_CONSTEXPR_AFTER_CXX14"},
{"lineNum":"   38","line":"typename iterator_traits<_RandIter>::difference_type"},
{"lineNum":"   39","line":"__distance(_RandIter __first, _RandIter __last, random_access_iterator_tag)"},
{"lineNum":"   40","line":"{"},
{"lineNum":"   41","line":"    return __last - __first;","class":"linePartCov","hits":"1","order":"660","possible_hits":"56",},
{"lineNum":"   42","line":"}"},
{"lineNum":"   43","line":""},
{"lineNum":"   44","line":"template <class _InputIter>"},
{"lineNum":"   45","line":"inline _LIBCPP_INLINE_VISIBILITY _LIBCPP_CONSTEXPR_AFTER_CXX14"},
{"lineNum":"   46","line":"typename iterator_traits<_InputIter>::difference_type"},
{"lineNum":"   47","line":"distance(_InputIter __first, _InputIter __last)"},
{"lineNum":"   48","line":"{"},
{"lineNum":"   49","line":"    return _VSTD::__distance(__first, __last, typename iterator_traits<_InputIter>::iterator_category());"},
{"lineNum":"   50","line":"}"},
{"lineNum":"   51","line":""},
{"lineNum":"   52","line":"_LIBCPP_END_NAMESPACE_STD"},
{"lineNum":"   53","line":""},
{"lineNum":"   54","line":"_LIBCPP_POP_MACROS"},
{"lineNum":"   55","line":""},
{"lineNum":"   56","line":"#endif // _LIBCPP___ITERATOR_DISTANCE_H"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 1, "covered" : 1,};
var merged_data = [];

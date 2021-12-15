var data = {lines:[
{"lineNum":"    1","line":"//===----------------------------------------------------------------------===//"},
{"lineNum":"    2","line":"//"},
{"lineNum":"    3","line":"// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions."},
{"lineNum":"    4","line":"// See https://llvm.org/LICENSE.txt for license information."},
{"lineNum":"    5","line":"// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception"},
{"lineNum":"    6","line":"//"},
{"lineNum":"    7","line":"//===----------------------------------------------------------------------===//"},
{"lineNum":"    8","line":""},
{"lineNum":"    9","line":"#ifndef _LIBCPP___ALGORITHM_TRANSFORM_H"},
{"lineNum":"   10","line":"#define _LIBCPP___ALGORITHM_TRANSFORM_H"},
{"lineNum":"   11","line":""},
{"lineNum":"   12","line":"#include <__config>"},
{"lineNum":"   13","line":""},
{"lineNum":"   14","line":"#if !defined(_LIBCPP_HAS_NO_PRAGMA_SYSTEM_HEADER)"},
{"lineNum":"   15","line":"#pragma GCC system_header"},
{"lineNum":"   16","line":"#endif"},
{"lineNum":"   17","line":""},
{"lineNum":"   18","line":"_LIBCPP_PUSH_MACROS"},
{"lineNum":"   19","line":"#include <__undef_macros>"},
{"lineNum":"   20","line":""},
{"lineNum":"   21","line":"_LIBCPP_BEGIN_NAMESPACE_STD"},
{"lineNum":"   22","line":""},
{"lineNum":"   23","line":"template <class _InputIterator, class _OutputIterator, class _UnaryOperation>"},
{"lineNum":"   24","line":"inline _LIBCPP_INLINE_VISIBILITY _LIBCPP_CONSTEXPR_AFTER_CXX17"},
{"lineNum":"   25","line":"_OutputIterator"},
{"lineNum":"   26","line":"transform(_InputIterator __first, _InputIterator __last, _OutputIterator __result, _UnaryOperation __op)"},
{"lineNum":"   27","line":"{"},
{"lineNum":"   28","line":"    for (; __first != __last; ++__first, (void) ++__result)","class":"lineNoCov","hits":"0","possible_hits":"435",},
{"lineNum":"   29","line":"        *__result = __op(*__first);","class":"lineNoCov","hits":"0","possible_hits":"770",},
{"lineNum":"   30","line":"    return __result;"},
{"lineNum":"   31","line":"}"},
{"lineNum":"   32","line":""},
{"lineNum":"   33","line":"template <class _InputIterator1, class _InputIterator2, class _OutputIterator, class _BinaryOperation>"},
{"lineNum":"   34","line":"inline _LIBCPP_INLINE_VISIBILITY _LIBCPP_CONSTEXPR_AFTER_CXX17"},
{"lineNum":"   35","line":"_OutputIterator"},
{"lineNum":"   36","line":"transform(_InputIterator1 __first1, _InputIterator1 __last1, _InputIterator2 __first2,"},
{"lineNum":"   37","line":"          _OutputIterator __result, _BinaryOperation __binary_op)"},
{"lineNum":"   38","line":"{"},
{"lineNum":"   39","line":"    for (; __first1 != __last1; ++__first1, (void) ++__first2, ++__result)"},
{"lineNum":"   40","line":"        *__result = __binary_op(*__first1, *__first2);"},
{"lineNum":"   41","line":"    return __result;"},
{"lineNum":"   42","line":"}"},
{"lineNum":"   43","line":""},
{"lineNum":"   44","line":"_LIBCPP_END_NAMESPACE_STD"},
{"lineNum":"   45","line":""},
{"lineNum":"   46","line":"_LIBCPP_POP_MACROS"},
{"lineNum":"   47","line":""},
{"lineNum":"   48","line":"#endif // _LIBCPP___ALGORITHM_TRANSFORM_H"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 2, "covered" : 0,};
var merged_data = [];

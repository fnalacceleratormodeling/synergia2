var data = {lines:[
{"lineNum":"    1","line":"//===----------------------------------------------------------------------===//"},
{"lineNum":"    2","line":"//"},
{"lineNum":"    3","line":"// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions."},
{"lineNum":"    4","line":"// See https://llvm.org/LICENSE.txt for license information."},
{"lineNum":"    5","line":"// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception"},
{"lineNum":"    6","line":"//"},
{"lineNum":"    7","line":"//===----------------------------------------------------------------------===//"},
{"lineNum":"    8","line":""},
{"lineNum":"    9","line":"#ifndef _LIBCPP___ALGORITHM_SORT_H"},
{"lineNum":"   10","line":"#define _LIBCPP___ALGORITHM_SORT_H"},
{"lineNum":"   11","line":""},
{"lineNum":"   12","line":"#include <__config>"},
{"lineNum":"   13","line":"#include <__algorithm/comp.h>"},
{"lineNum":"   14","line":"#include <__algorithm/comp_ref_type.h>"},
{"lineNum":"   15","line":"#include <__algorithm/min_element.h>"},
{"lineNum":"   16","line":"#include <__algorithm/partial_sort.h>"},
{"lineNum":"   17","line":"#include <__algorithm/unwrap_iter.h>"},
{"lineNum":"   18","line":"#include <__utility/swap.h>"},
{"lineNum":"   19","line":"#include <memory>"},
{"lineNum":"   20","line":"#include <type_traits> // swap"},
{"lineNum":"   21","line":""},
{"lineNum":"   22","line":"#if !defined(_LIBCPP_HAS_NO_PRAGMA_SYSTEM_HEADER)"},
{"lineNum":"   23","line":"#pragma GCC system_header"},
{"lineNum":"   24","line":"#endif"},
{"lineNum":"   25","line":""},
{"lineNum":"   26","line":"_LIBCPP_PUSH_MACROS"},
{"lineNum":"   27","line":"#include <__undef_macros>"},
{"lineNum":"   28","line":""},
{"lineNum":"   29","line":"_LIBCPP_BEGIN_NAMESPACE_STD"},
{"lineNum":"   30","line":""},
{"lineNum":"   31","line":"// stable, 2-3 compares, 0-2 swaps"},
{"lineNum":"   32","line":""},
{"lineNum":"   33","line":"template <class _Compare, class _ForwardIterator>"},
{"lineNum":"   34","line":"_LIBCPP_CONSTEXPR_AFTER_CXX11 unsigned"},
{"lineNum":"   35","line":"__sort3(_ForwardIterator __x, _ForwardIterator __y, _ForwardIterator __z, _Compare __c)"},
{"lineNum":"   36","line":"{"},
{"lineNum":"   37","line":"    unsigned __r = 0;"},
{"lineNum":"   38","line":"    if (!__c(*__y, *__x))          // if x <= y"},
{"lineNum":"   39","line":"    {"},
{"lineNum":"   40","line":"        if (!__c(*__z, *__y))      // if y <= z"},
{"lineNum":"   41","line":"            return __r;            // x <= y && y <= z"},
{"lineNum":"   42","line":"                                   // x <= y && y > z"},
{"lineNum":"   43","line":"        swap(*__y, *__z);          // x <= z && y < z"},
{"lineNum":"   44","line":"        __r = 1;"},
{"lineNum":"   45","line":"        if (__c(*__y, *__x))       // if x > y"},
{"lineNum":"   46","line":"        {"},
{"lineNum":"   47","line":"            swap(*__x, *__y);      // x < y && y <= z"},
{"lineNum":"   48","line":"            __r = 2;"},
{"lineNum":"   49","line":"        }"},
{"lineNum":"   50","line":"        return __r;                // x <= y && y < z"},
{"lineNum":"   51","line":"    }"},
{"lineNum":"   52","line":"    if (__c(*__z, *__y))           // x > y, if y > z"},
{"lineNum":"   53","line":"    {"},
{"lineNum":"   54","line":"        swap(*__x, *__z);          // x < y && y < z"},
{"lineNum":"   55","line":"        __r = 1;"},
{"lineNum":"   56","line":"        return __r;"},
{"lineNum":"   57","line":"    }"},
{"lineNum":"   58","line":"    swap(*__x, *__y);              // x > y && y <= z"},
{"lineNum":"   59","line":"    __r = 1;                       // x < y && x <= z"},
{"lineNum":"   60","line":"    if (__c(*__z, *__y))           // if y > z"},
{"lineNum":"   61","line":"    {"},
{"lineNum":"   62","line":"        swap(*__y, *__z);          // x <= y && y < z"},
{"lineNum":"   63","line":"        __r = 2;"},
{"lineNum":"   64","line":"    }"},
{"lineNum":"   65","line":"    return __r;"},
{"lineNum":"   66","line":"}                                  // x <= y && y <= z"},
{"lineNum":"   67","line":""},
{"lineNum":"   68","line":"// stable, 3-6 compares, 0-5 swaps"},
{"lineNum":"   69","line":""},
{"lineNum":"   70","line":"template <class _Compare, class _ForwardIterator>"},
{"lineNum":"   71","line":"unsigned"},
{"lineNum":"   72","line":"__sort4(_ForwardIterator __x1, _ForwardIterator __x2, _ForwardIterator __x3,"},
{"lineNum":"   73","line":"            _ForwardIterator __x4, _Compare __c)"},
{"lineNum":"   74","line":"{"},
{"lineNum":"   75","line":"    unsigned __r = _VSTD::__sort3<_Compare>(__x1, __x2, __x3, __c);"},
{"lineNum":"   76","line":"    if (__c(*__x4, *__x3))"},
{"lineNum":"   77","line":"    {"},
{"lineNum":"   78","line":"        swap(*__x3, *__x4);"},
{"lineNum":"   79","line":"        ++__r;"},
{"lineNum":"   80","line":"        if (__c(*__x3, *__x2))"},
{"lineNum":"   81","line":"        {"},
{"lineNum":"   82","line":"            swap(*__x2, *__x3);"},
{"lineNum":"   83","line":"            ++__r;"},
{"lineNum":"   84","line":"            if (__c(*__x2, *__x1))"},
{"lineNum":"   85","line":"            {"},
{"lineNum":"   86","line":"                swap(*__x1, *__x2);"},
{"lineNum":"   87","line":"                ++__r;"},
{"lineNum":"   88","line":"            }"},
{"lineNum":"   89","line":"        }"},
{"lineNum":"   90","line":"    }"},
{"lineNum":"   91","line":"    return __r;"},
{"lineNum":"   92","line":"}"},
{"lineNum":"   93","line":""},
{"lineNum":"   94","line":"// stable, 4-10 compares, 0-9 swaps"},
{"lineNum":"   95","line":""},
{"lineNum":"   96","line":"template <class _Compare, class _ForwardIterator>"},
{"lineNum":"   97","line":"_LIBCPP_HIDDEN"},
{"lineNum":"   98","line":"unsigned"},
{"lineNum":"   99","line":"__sort5(_ForwardIterator __x1, _ForwardIterator __x2, _ForwardIterator __x3,"},
{"lineNum":"  100","line":"            _ForwardIterator __x4, _ForwardIterator __x5, _Compare __c)"},
{"lineNum":"  101","line":"{"},
{"lineNum":"  102","line":"    unsigned __r = _VSTD::__sort4<_Compare>(__x1, __x2, __x3, __x4, __c);"},
{"lineNum":"  103","line":"    if (__c(*__x5, *__x4))"},
{"lineNum":"  104","line":"    {"},
{"lineNum":"  105","line":"        swap(*__x4, *__x5);"},
{"lineNum":"  106","line":"        ++__r;"},
{"lineNum":"  107","line":"        if (__c(*__x4, *__x3))"},
{"lineNum":"  108","line":"        {"},
{"lineNum":"  109","line":"            swap(*__x3, *__x4);"},
{"lineNum":"  110","line":"            ++__r;"},
{"lineNum":"  111","line":"            if (__c(*__x3, *__x2))"},
{"lineNum":"  112","line":"            {"},
{"lineNum":"  113","line":"                swap(*__x2, *__x3);"},
{"lineNum":"  114","line":"                ++__r;"},
{"lineNum":"  115","line":"                if (__c(*__x2, *__x1))"},
{"lineNum":"  116","line":"                {"},
{"lineNum":"  117","line":"                    swap(*__x1, *__x2);"},
{"lineNum":"  118","line":"                    ++__r;"},
{"lineNum":"  119","line":"                }"},
{"lineNum":"  120","line":"            }"},
{"lineNum":"  121","line":"        }"},
{"lineNum":"  122","line":"    }"},
{"lineNum":"  123","line":"    return __r;"},
{"lineNum":"  124","line":"}"},
{"lineNum":"  125","line":""},
{"lineNum":"  126","line":"// Assumes size > 0"},
{"lineNum":"  127","line":"template <class _Compare, class _BidirectionalIterator>"},
{"lineNum":"  128","line":"_LIBCPP_CONSTEXPR_AFTER_CXX11 void"},
{"lineNum":"  129","line":"__selection_sort(_BidirectionalIterator __first, _BidirectionalIterator __last, _Compare __comp)"},
{"lineNum":"  130","line":"{"},
{"lineNum":"  131","line":"    _BidirectionalIterator __lm1 = __last;"},
{"lineNum":"  132","line":"    for (--__lm1; __first != __lm1; ++__first)"},
{"lineNum":"  133","line":"    {"},
{"lineNum":"  134","line":"        _BidirectionalIterator __i = _VSTD::min_element<_BidirectionalIterator,"},
{"lineNum":"  135","line":"                                                        typename add_lvalue_reference<_Compare>::type>"},
{"lineNum":"  136","line":"                                                       (__first, __last, __comp);"},
{"lineNum":"  137","line":"        if (__i != __first)"},
{"lineNum":"  138","line":"            swap(*__first, *__i);"},
{"lineNum":"  139","line":"    }"},
{"lineNum":"  140","line":"}"},
{"lineNum":"  141","line":""},
{"lineNum":"  142","line":"template <class _Compare, class _BidirectionalIterator>"},
{"lineNum":"  143","line":"void"},
{"lineNum":"  144","line":"__insertion_sort(_BidirectionalIterator __first, _BidirectionalIterator __last, _Compare __comp)"},
{"lineNum":"  145","line":"{"},
{"lineNum":"  146","line":"    typedef typename iterator_traits<_BidirectionalIterator>::value_type value_type;"},
{"lineNum":"  147","line":"    if (__first != __last)"},
{"lineNum":"  148","line":"    {"},
{"lineNum":"  149","line":"        _BidirectionalIterator __i = __first;"},
{"lineNum":"  150","line":"        for (++__i; __i != __last; ++__i)"},
{"lineNum":"  151","line":"        {"},
{"lineNum":"  152","line":"            _BidirectionalIterator __j = __i;"},
{"lineNum":"  153","line":"            value_type __t(_VSTD::move(*__j));"},
{"lineNum":"  154","line":"            for (_BidirectionalIterator __k = __i; __k != __first && __comp(__t,  *--__k); --__j)"},
{"lineNum":"  155","line":"                *__j = _VSTD::move(*__k);"},
{"lineNum":"  156","line":"            *__j = _VSTD::move(__t);"},
{"lineNum":"  157","line":"        }"},
{"lineNum":"  158","line":"    }"},
{"lineNum":"  159","line":"}"},
{"lineNum":"  160","line":""},
{"lineNum":"  161","line":"template <class _Compare, class _RandomAccessIterator>"},
{"lineNum":"  162","line":"void"},
{"lineNum":"  163","line":"__insertion_sort_3(_RandomAccessIterator __first, _RandomAccessIterator __last, _Compare __comp)"},
{"lineNum":"  164","line":"{"},
{"lineNum":"  165","line":"    typedef typename iterator_traits<_RandomAccessIterator>::value_type value_type;"},
{"lineNum":"  166","line":"    _RandomAccessIterator __j = __first+2;"},
{"lineNum":"  167","line":"    _VSTD::__sort3<_Compare>(__first, __first+1, __j, __comp);"},
{"lineNum":"  168","line":"    for (_RandomAccessIterator __i = __j+1; __i != __last; ++__i)"},
{"lineNum":"  169","line":"    {"},
{"lineNum":"  170","line":"        if (__comp(*__i, *__j))"},
{"lineNum":"  171","line":"        {"},
{"lineNum":"  172","line":"            value_type __t(_VSTD::move(*__i));"},
{"lineNum":"  173","line":"            _RandomAccessIterator __k = __j;"},
{"lineNum":"  174","line":"            __j = __i;"},
{"lineNum":"  175","line":"            do"},
{"lineNum":"  176","line":"            {"},
{"lineNum":"  177","line":"                *__j = _VSTD::move(*__k);"},
{"lineNum":"  178","line":"                __j = __k;"},
{"lineNum":"  179","line":"            } while (__j != __first && __comp(__t, *--__k));"},
{"lineNum":"  180","line":"            *__j = _VSTD::move(__t);"},
{"lineNum":"  181","line":"        }"},
{"lineNum":"  182","line":"        __j = __i;"},
{"lineNum":"  183","line":"    }"},
{"lineNum":"  184","line":"}"},
{"lineNum":"  185","line":""},
{"lineNum":"  186","line":"template <class _Compare, class _RandomAccessIterator>"},
{"lineNum":"  187","line":"bool"},
{"lineNum":"  188","line":"__insertion_sort_incomplete(_RandomAccessIterator __first, _RandomAccessIterator __last, _Compare __comp)"},
{"lineNum":"  189","line":"{"},
{"lineNum":"  190","line":"    switch (__last - __first)"},
{"lineNum":"  191","line":"    {"},
{"lineNum":"  192","line":"    case 0:"},
{"lineNum":"  193","line":"    case 1:"},
{"lineNum":"  194","line":"        return true;"},
{"lineNum":"  195","line":"    case 2:"},
{"lineNum":"  196","line":"        if (__comp(*--__last, *__first))"},
{"lineNum":"  197","line":"            swap(*__first, *__last);"},
{"lineNum":"  198","line":"        return true;"},
{"lineNum":"  199","line":"    case 3:"},
{"lineNum":"  200","line":"        _VSTD::__sort3<_Compare>(__first, __first+1, --__last, __comp);"},
{"lineNum":"  201","line":"        return true;"},
{"lineNum":"  202","line":"    case 4:"},
{"lineNum":"  203","line":"        _VSTD::__sort4<_Compare>(__first, __first+1, __first+2, --__last, __comp);"},
{"lineNum":"  204","line":"        return true;"},
{"lineNum":"  205","line":"    case 5:"},
{"lineNum":"  206","line":"        _VSTD::__sort5<_Compare>(__first, __first+1, __first+2, __first+3, --__last, __comp);"},
{"lineNum":"  207","line":"        return true;"},
{"lineNum":"  208","line":"    }"},
{"lineNum":"  209","line":"    typedef typename iterator_traits<_RandomAccessIterator>::value_type value_type;"},
{"lineNum":"  210","line":"    _RandomAccessIterator __j = __first+2;"},
{"lineNum":"  211","line":"    _VSTD::__sort3<_Compare>(__first, __first+1, __j, __comp);"},
{"lineNum":"  212","line":"    const unsigned __limit = 8;"},
{"lineNum":"  213","line":"    unsigned __count = 0;"},
{"lineNum":"  214","line":"    for (_RandomAccessIterator __i = __j+1; __i != __last; ++__i)"},
{"lineNum":"  215","line":"    {"},
{"lineNum":"  216","line":"        if (__comp(*__i, *__j))"},
{"lineNum":"  217","line":"        {"},
{"lineNum":"  218","line":"            value_type __t(_VSTD::move(*__i));"},
{"lineNum":"  219","line":"            _RandomAccessIterator __k = __j;"},
{"lineNum":"  220","line":"            __j = __i;"},
{"lineNum":"  221","line":"            do"},
{"lineNum":"  222","line":"            {"},
{"lineNum":"  223","line":"                *__j = _VSTD::move(*__k);"},
{"lineNum":"  224","line":"                __j = __k;"},
{"lineNum":"  225","line":"            } while (__j != __first && __comp(__t, *--__k));"},
{"lineNum":"  226","line":"            *__j = _VSTD::move(__t);"},
{"lineNum":"  227","line":"            if (++__count == __limit)"},
{"lineNum":"  228","line":"                return ++__i == __last;"},
{"lineNum":"  229","line":"        }"},
{"lineNum":"  230","line":"        __j = __i;"},
{"lineNum":"  231","line":"    }"},
{"lineNum":"  232","line":"    return true;"},
{"lineNum":"  233","line":"}"},
{"lineNum":"  234","line":""},
{"lineNum":"  235","line":"template <class _Compare, class _BidirectionalIterator>"},
{"lineNum":"  236","line":"void"},
{"lineNum":"  237","line":"__insertion_sort_move(_BidirectionalIterator __first1, _BidirectionalIterator __last1,"},
{"lineNum":"  238","line":"                      typename iterator_traits<_BidirectionalIterator>::value_type* __first2, _Compare __comp)"},
{"lineNum":"  239","line":"{"},
{"lineNum":"  240","line":"    typedef typename iterator_traits<_BidirectionalIterator>::value_type value_type;"},
{"lineNum":"  241","line":"    if (__first1 != __last1)"},
{"lineNum":"  242","line":"    {"},
{"lineNum":"  243","line":"        __destruct_n __d(0);"},
{"lineNum":"  244","line":"        unique_ptr<value_type, __destruct_n&> __h(__first2, __d);"},
{"lineNum":"  245","line":"        value_type* __last2 = __first2;"},
{"lineNum":"  246","line":"        ::new ((void*)__last2) value_type(_VSTD::move(*__first1));"},
{"lineNum":"  247","line":"        __d.template __incr<value_type>();"},
{"lineNum":"  248","line":"        for (++__last2; ++__first1 != __last1; ++__last2)"},
{"lineNum":"  249","line":"        {"},
{"lineNum":"  250","line":"            value_type* __j2 = __last2;"},
{"lineNum":"  251","line":"            value_type* __i2 = __j2;"},
{"lineNum":"  252","line":"            if (__comp(*__first1, *--__i2))"},
{"lineNum":"  253","line":"            {"},
{"lineNum":"  254","line":"                ::new ((void*)__j2) value_type(_VSTD::move(*__i2));"},
{"lineNum":"  255","line":"                __d.template __incr<value_type>();"},
{"lineNum":"  256","line":"                for (--__j2; __i2 != __first2 && __comp(*__first1,  *--__i2); --__j2)"},
{"lineNum":"  257","line":"                    *__j2 = _VSTD::move(*__i2);"},
{"lineNum":"  258","line":"                *__j2 = _VSTD::move(*__first1);"},
{"lineNum":"  259","line":"            }"},
{"lineNum":"  260","line":"            else"},
{"lineNum":"  261","line":"            {"},
{"lineNum":"  262","line":"                ::new ((void*)__j2) value_type(_VSTD::move(*__first1));"},
{"lineNum":"  263","line":"                __d.template __incr<value_type>();"},
{"lineNum":"  264","line":"            }"},
{"lineNum":"  265","line":"        }"},
{"lineNum":"  266","line":"        __h.release();"},
{"lineNum":"  267","line":"    }"},
{"lineNum":"  268","line":"}"},
{"lineNum":"  269","line":""},
{"lineNum":"  270","line":"template <class _Compare, class _RandomAccessIterator>"},
{"lineNum":"  271","line":"void"},
{"lineNum":"  272","line":"__sort(_RandomAccessIterator __first, _RandomAccessIterator __last, _Compare __comp)"},
{"lineNum":"  273","line":"{"},
{"lineNum":"  274","line":"    typedef typename iterator_traits<_RandomAccessIterator>::difference_type difference_type;"},
{"lineNum":"  275","line":"    typedef typename iterator_traits<_RandomAccessIterator>::value_type value_type;"},
{"lineNum":"  276","line":"    const difference_type __limit = is_trivially_copy_constructible<value_type>::value &&"},
{"lineNum":"  277","line":"                                    is_trivially_copy_assignable<value_type>::value ? 30 : 6;"},
{"lineNum":"  278","line":"    while (true)"},
{"lineNum":"  279","line":"    {"},
{"lineNum":"  280","line":"    __restart:"},
{"lineNum":"  281","line":"        difference_type __len = __last - __first;"},
{"lineNum":"  282","line":"        switch (__len)"},
{"lineNum":"  283","line":"        {"},
{"lineNum":"  284","line":"        case 0:"},
{"lineNum":"  285","line":"        case 1:"},
{"lineNum":"  286","line":"            return;"},
{"lineNum":"  287","line":"        case 2:"},
{"lineNum":"  288","line":"            if (__comp(*--__last, *__first))"},
{"lineNum":"  289","line":"                swap(*__first, *__last);"},
{"lineNum":"  290","line":"            return;"},
{"lineNum":"  291","line":"        case 3:"},
{"lineNum":"  292","line":"            _VSTD::__sort3<_Compare>(__first, __first+1, --__last, __comp);"},
{"lineNum":"  293","line":"            return;"},
{"lineNum":"  294","line":"        case 4:"},
{"lineNum":"  295","line":"            _VSTD::__sort4<_Compare>(__first, __first+1, __first+2, --__last, __comp);"},
{"lineNum":"  296","line":"            return;"},
{"lineNum":"  297","line":"        case 5:"},
{"lineNum":"  298","line":"            _VSTD::__sort5<_Compare>(__first, __first+1, __first+2, __first+3, --__last, __comp);"},
{"lineNum":"  299","line":"            return;"},
{"lineNum":"  300","line":"        }"},
{"lineNum":"  301","line":"        if (__len <= __limit)"},
{"lineNum":"  302","line":"        {"},
{"lineNum":"  303","line":"            _VSTD::__insertion_sort_3<_Compare>(__first, __last, __comp);"},
{"lineNum":"  304","line":"            return;"},
{"lineNum":"  305","line":"        }"},
{"lineNum":"  306","line":"        // __len > 5"},
{"lineNum":"  307","line":"        _RandomAccessIterator __m = __first;"},
{"lineNum":"  308","line":"        _RandomAccessIterator __lm1 = __last;"},
{"lineNum":"  309","line":"        --__lm1;"},
{"lineNum":"  310","line":"        unsigned __n_swaps;"},
{"lineNum":"  311","line":"        {"},
{"lineNum":"  312","line":"        difference_type __delta;"},
{"lineNum":"  313","line":"        if (__len >= 1000)"},
{"lineNum":"  314","line":"        {"},
{"lineNum":"  315","line":"            __delta = __len/2;"},
{"lineNum":"  316","line":"            __m += __delta;"},
{"lineNum":"  317","line":"            __delta /= 2;"},
{"lineNum":"  318","line":"            __n_swaps = _VSTD::__sort5<_Compare>(__first, __first + __delta, __m, __m+__delta, __lm1, __comp);"},
{"lineNum":"  319","line":"        }"},
{"lineNum":"  320","line":"        else"},
{"lineNum":"  321","line":"        {"},
{"lineNum":"  322","line":"            __delta = __len/2;"},
{"lineNum":"  323","line":"            __m += __delta;"},
{"lineNum":"  324","line":"            __n_swaps = _VSTD::__sort3<_Compare>(__first, __m, __lm1, __comp);"},
{"lineNum":"  325","line":"        }"},
{"lineNum":"  326","line":"        }"},
{"lineNum":"  327","line":"        // *__m is median"},
{"lineNum":"  328","line":"        // partition [__first, __m) < *__m and *__m <= [__m, __last)"},
{"lineNum":"  329","line":"        // (this inhibits tossing elements equivalent to __m around unnecessarily)"},
{"lineNum":"  330","line":"        _RandomAccessIterator __i = __first;"},
{"lineNum":"  331","line":"        _RandomAccessIterator __j = __lm1;"},
{"lineNum":"  332","line":"        // j points beyond range to be tested, *__m is known to be <= *__lm1"},
{"lineNum":"  333","line":"        // The search going up is known to be guarded but the search coming down isn\'t."},
{"lineNum":"  334","line":"        // Prime the downward search with a guard."},
{"lineNum":"  335","line":"        if (!__comp(*__i, *__m))  // if *__first == *__m"},
{"lineNum":"  336","line":"        {"},
{"lineNum":"  337","line":"            // *__first == *__m, *__first doesn\'t go in first part"},
{"lineNum":"  338","line":"            // manually guard downward moving __j against __i"},
{"lineNum":"  339","line":"            while (true)"},
{"lineNum":"  340","line":"            {"},
{"lineNum":"  341","line":"                if (__i == --__j)"},
{"lineNum":"  342","line":"                {"},
{"lineNum":"  343","line":"                    // *__first == *__m, *__m <= all other elements"},
{"lineNum":"  344","line":"                    // Parition instead into [__first, __i) == *__first and *__first < [__i, __last)"},
{"lineNum":"  345","line":"                    ++__i;  // __first + 1"},
{"lineNum":"  346","line":"                    __j = __last;"},
{"lineNum":"  347","line":"                    if (!__comp(*__first, *--__j))  // we need a guard if *__first == *(__last-1)"},
{"lineNum":"  348","line":"                    {"},
{"lineNum":"  349","line":"                        while (true)"},
{"lineNum":"  350","line":"                        {"},
{"lineNum":"  351","line":"                            if (__i == __j)"},
{"lineNum":"  352","line":"                                return;  // [__first, __last) all equivalent elements"},
{"lineNum":"  353","line":"                            if (__comp(*__first, *__i))"},
{"lineNum":"  354","line":"                            {"},
{"lineNum":"  355","line":"                                swap(*__i, *__j);"},
{"lineNum":"  356","line":"                                ++__n_swaps;"},
{"lineNum":"  357","line":"                                ++__i;"},
{"lineNum":"  358","line":"                                break;"},
{"lineNum":"  359","line":"                            }"},
{"lineNum":"  360","line":"                            ++__i;"},
{"lineNum":"  361","line":"                        }"},
{"lineNum":"  362","line":"                    }"},
{"lineNum":"  363","line":"                    // [__first, __i) == *__first and *__first < [__j, __last) and __j == __last - 1"},
{"lineNum":"  364","line":"                    if (__i == __j)"},
{"lineNum":"  365","line":"                        return;"},
{"lineNum":"  366","line":"                    while (true)"},
{"lineNum":"  367","line":"                    {"},
{"lineNum":"  368","line":"                        while (!__comp(*__first, *__i))"},
{"lineNum":"  369","line":"                            ++__i;"},
{"lineNum":"  370","line":"                        while (__comp(*__first, *--__j))"},
{"lineNum":"  371","line":"                            ;"},
{"lineNum":"  372","line":"                        if (__i >= __j)"},
{"lineNum":"  373","line":"                            break;"},
{"lineNum":"  374","line":"                        swap(*__i, *__j);"},
{"lineNum":"  375","line":"                        ++__n_swaps;"},
{"lineNum":"  376","line":"                        ++__i;"},
{"lineNum":"  377","line":"                    }"},
{"lineNum":"  378","line":"                    // [__first, __i) == *__first and *__first < [__i, __last)"},
{"lineNum":"  379","line":"                    // The first part is sorted, sort the second part"},
{"lineNum":"  380","line":"                    // _VSTD::__sort<_Compare>(__i, __last, __comp);"},
{"lineNum":"  381","line":"                    __first = __i;"},
{"lineNum":"  382","line":"                    goto __restart;"},
{"lineNum":"  383","line":"                }"},
{"lineNum":"  384","line":"                if (__comp(*__j, *__m))"},
{"lineNum":"  385","line":"                {"},
{"lineNum":"  386","line":"                    swap(*__i, *__j);"},
{"lineNum":"  387","line":"                    ++__n_swaps;"},
{"lineNum":"  388","line":"                    break;  // found guard for downward moving __j, now use unguarded partition"},
{"lineNum":"  389","line":"                }"},
{"lineNum":"  390","line":"            }"},
{"lineNum":"  391","line":"        }"},
{"lineNum":"  392","line":"        // It is known that *__i < *__m"},
{"lineNum":"  393","line":"        ++__i;"},
{"lineNum":"  394","line":"        // j points beyond range to be tested, *__m is known to be <= *__lm1"},
{"lineNum":"  395","line":"        // if not yet partitioned..."},
{"lineNum":"  396","line":"        if (__i < __j)"},
{"lineNum":"  397","line":"        {"},
{"lineNum":"  398","line":"            // known that *(__i - 1) < *__m"},
{"lineNum":"  399","line":"            // known that __i <= __m"},
{"lineNum":"  400","line":"            while (true)"},
{"lineNum":"  401","line":"            {"},
{"lineNum":"  402","line":"                // __m still guards upward moving __i"},
{"lineNum":"  403","line":"                while (__comp(*__i, *__m))"},
{"lineNum":"  404","line":"                    ++__i;"},
{"lineNum":"  405","line":"                // It is now known that a guard exists for downward moving __j"},
{"lineNum":"  406","line":"                while (!__comp(*--__j, *__m))"},
{"lineNum":"  407","line":"                    ;"},
{"lineNum":"  408","line":"                if (__i > __j)"},
{"lineNum":"  409","line":"                    break;"},
{"lineNum":"  410","line":"                swap(*__i, *__j);"},
{"lineNum":"  411","line":"                ++__n_swaps;"},
{"lineNum":"  412","line":"                // It is known that __m != __j"},
{"lineNum":"  413","line":"                // If __m just moved, follow it"},
{"lineNum":"  414","line":"                if (__m == __i)"},
{"lineNum":"  415","line":"                    __m = __j;"},
{"lineNum":"  416","line":"                ++__i;"},
{"lineNum":"  417","line":"            }"},
{"lineNum":"  418","line":"        }"},
{"lineNum":"  419","line":"        // [__first, __i) < *__m and *__m <= [__i, __last)"},
{"lineNum":"  420","line":"        if (__i != __m && __comp(*__m, *__i))"},
{"lineNum":"  421","line":"        {"},
{"lineNum":"  422","line":"            swap(*__i, *__m);"},
{"lineNum":"  423","line":"            ++__n_swaps;"},
{"lineNum":"  424","line":"        }"},
{"lineNum":"  425","line":"        // [__first, __i) < *__i and *__i <= [__i+1, __last)"},
{"lineNum":"  426","line":"        // If we were given a perfect partition, see if insertion sort is quick..."},
{"lineNum":"  427","line":"        if (__n_swaps == 0)"},
{"lineNum":"  428","line":"        {"},
{"lineNum":"  429","line":"            bool __fs = _VSTD::__insertion_sort_incomplete<_Compare>(__first, __i, __comp);"},
{"lineNum":"  430","line":"            if (_VSTD::__insertion_sort_incomplete<_Compare>(__i+1, __last, __comp))"},
{"lineNum":"  431","line":"            {"},
{"lineNum":"  432","line":"                if (__fs)"},
{"lineNum":"  433","line":"                    return;"},
{"lineNum":"  434","line":"                __last = __i;"},
{"lineNum":"  435","line":"                continue;"},
{"lineNum":"  436","line":"            }"},
{"lineNum":"  437","line":"            else"},
{"lineNum":"  438","line":"            {"},
{"lineNum":"  439","line":"                if (__fs)"},
{"lineNum":"  440","line":"                {"},
{"lineNum":"  441","line":"                    __first = ++__i;"},
{"lineNum":"  442","line":"                    continue;"},
{"lineNum":"  443","line":"                }"},
{"lineNum":"  444","line":"            }"},
{"lineNum":"  445","line":"        }"},
{"lineNum":"  446","line":"        // sort smaller range with recursive call and larger with tail recursion elimination"},
{"lineNum":"  447","line":"        if (__i - __first < __last - __i)"},
{"lineNum":"  448","line":"        {"},
{"lineNum":"  449","line":"            _VSTD::__sort<_Compare>(__first, __i, __comp);"},
{"lineNum":"  450","line":"            // _VSTD::__sort<_Compare>(__i+1, __last, __comp);"},
{"lineNum":"  451","line":"            __first = ++__i;"},
{"lineNum":"  452","line":"        }"},
{"lineNum":"  453","line":"        else"},
{"lineNum":"  454","line":"        {"},
{"lineNum":"  455","line":"            _VSTD::__sort<_Compare>(__i+1, __last, __comp);"},
{"lineNum":"  456","line":"            // _VSTD::__sort<_Compare>(__first, __i, __comp);"},
{"lineNum":"  457","line":"            __last = __i;"},
{"lineNum":"  458","line":"        }"},
{"lineNum":"  459","line":"    }"},
{"lineNum":"  460","line":"}"},
{"lineNum":"  461","line":""},
{"lineNum":"  462","line":"template <class _Compare, class _Tp>"},
{"lineNum":"  463","line":"inline _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  464","line":"void"},
{"lineNum":"  465","line":"__sort(_Tp** __first, _Tp** __last, __less<_Tp*>&)"},
{"lineNum":"  466","line":"{"},
{"lineNum":"  467","line":"    __less<uintptr_t> __comp;"},
{"lineNum":"  468","line":"    _VSTD::__sort<__less<uintptr_t>&, uintptr_t*>((uintptr_t*)__first, (uintptr_t*)__last, __comp);"},
{"lineNum":"  469","line":"}"},
{"lineNum":"  470","line":""},
{"lineNum":"  471","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<char>&, char*>(char*, char*, __less<char>&))"},
{"lineNum":"  472","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<wchar_t>&, wchar_t*>(wchar_t*, wchar_t*, __less<wchar_t>&))"},
{"lineNum":"  473","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<signed char>&, signed char*>(signed char*, signed char*, __less<signed char>&))"},
{"lineNum":"  474","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<unsigned char>&, unsigned char*>(unsigned char*, unsigned char*, __less<unsigned char>&))"},
{"lineNum":"  475","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<short>&, short*>(short*, short*, __less<short>&))"},
{"lineNum":"  476","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<unsigned short>&, unsigned short*>(unsigned short*, unsigned short*, __less<unsigned short>&))"},
{"lineNum":"  477","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<int>&, int*>(int*, int*, __less<int>&))"},
{"lineNum":"  478","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<unsigned>&, unsigned*>(unsigned*, unsigned*, __less<unsigned>&))"},
{"lineNum":"  479","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<long>&, long*>(long*, long*, __less<long>&))"},
{"lineNum":"  480","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<unsigned long>&, unsigned long*>(unsigned long*, unsigned long*, __less<unsigned long>&))"},
{"lineNum":"  481","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<long long>&, long long*>(long long*, long long*, __less<long long>&))"},
{"lineNum":"  482","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<unsigned long long>&, unsigned long long*>(unsigned long long*, unsigned long long*, __less<unsigned long long>&))"},
{"lineNum":"  483","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<float>&, float*>(float*, float*, __less<float>&))"},
{"lineNum":"  484","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<double>&, double*>(double*, double*, __less<double>&))"},
{"lineNum":"  485","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS void __sort<__less<long double>&, long double*>(long double*, long double*, __less<long double>&))"},
{"lineNum":"  486","line":""},
{"lineNum":"  487","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<char>&, char*>(char*, char*, __less<char>&))"},
{"lineNum":"  488","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<wchar_t>&, wchar_t*>(wchar_t*, wchar_t*, __less<wchar_t>&))"},
{"lineNum":"  489","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<signed char>&, signed char*>(signed char*, signed char*, __less<signed char>&))"},
{"lineNum":"  490","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<unsigned char>&, unsigned char*>(unsigned char*, unsigned char*, __less<unsigned char>&))"},
{"lineNum":"  491","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<short>&, short*>(short*, short*, __less<short>&))"},
{"lineNum":"  492","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<unsigned short>&, unsigned short*>(unsigned short*, unsigned short*, __less<unsigned short>&))"},
{"lineNum":"  493","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<int>&, int*>(int*, int*, __less<int>&))"},
{"lineNum":"  494","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<unsigned>&, unsigned*>(unsigned*, unsigned*, __less<unsigned>&))"},
{"lineNum":"  495","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<long>&, long*>(long*, long*, __less<long>&))"},
{"lineNum":"  496","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<unsigned long>&, unsigned long*>(unsigned long*, unsigned long*, __less<unsigned long>&))"},
{"lineNum":"  497","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<long long>&, long long*>(long long*, long long*, __less<long long>&))"},
{"lineNum":"  498","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<unsigned long long>&, unsigned long long*>(unsigned long long*, unsigned long long*, __less<unsigned long long>&))"},
{"lineNum":"  499","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<float>&, float*>(float*, float*, __less<float>&))"},
{"lineNum":"  500","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<double>&, double*>(double*, double*, __less<double>&))"},
{"lineNum":"  501","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS bool __insertion_sort_incomplete<__less<long double>&, long double*>(long double*, long double*, __less<long double>&))"},
{"lineNum":"  502","line":""},
{"lineNum":"  503","line":"_LIBCPP_EXTERN_TEMPLATE(_LIBCPP_FUNC_VIS unsigned __sort5<__less<long double>&, long double*>(long double*, long double*, long double*, long double*, long double*, __less<long double>&))"},
{"lineNum":"  504","line":""},
{"lineNum":"  505","line":"template <class _RandomAccessIterator, class _Compare>"},
{"lineNum":"  506","line":"inline _LIBCPP_INLINE_VISIBILITY _LIBCPP_CONSTEXPR_AFTER_CXX17"},
{"lineNum":"  507","line":"void"},
{"lineNum":"  508","line":"sort(_RandomAccessIterator __first, _RandomAccessIterator __last, _Compare __comp)"},
{"lineNum":"  509","line":"{"},
{"lineNum":"  510","line":"    typedef typename __comp_ref_type<_Compare>::type _Comp_ref;"},
{"lineNum":"  511","line":"    if (__libcpp_is_constant_evaluated()) {"},
{"lineNum":"  512","line":"        _VSTD::__partial_sort<_Comp_ref>(__first, __last, __last, _Comp_ref(__comp));"},
{"lineNum":"  513","line":"    } else {"},
{"lineNum":"  514","line":"        _VSTD::__sort<_Comp_ref>(_VSTD::__unwrap_iter(__first), _VSTD::__unwrap_iter(__last), _Comp_ref(__comp));","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  515","line":"    }"},
{"lineNum":"  516","line":"}"},
{"lineNum":"  517","line":""},
{"lineNum":"  518","line":"template <class _RandomAccessIterator>"},
{"lineNum":"  519","line":"inline _LIBCPP_INLINE_VISIBILITY _LIBCPP_CONSTEXPR_AFTER_CXX17"},
{"lineNum":"  520","line":"void"},
{"lineNum":"  521","line":"sort(_RandomAccessIterator __first, _RandomAccessIterator __last)"},
{"lineNum":"  522","line":"{"},
{"lineNum":"  523","line":"    _VSTD::sort(__first, __last, __less<typename iterator_traits<_RandomAccessIterator>::value_type>());"},
{"lineNum":"  524","line":"}"},
{"lineNum":"  525","line":""},
{"lineNum":"  526","line":"_LIBCPP_END_NAMESPACE_STD"},
{"lineNum":"  527","line":""},
{"lineNum":"  528","line":"_LIBCPP_POP_MACROS"},
{"lineNum":"  529","line":""},
{"lineNum":"  530","line":"#endif // _LIBCPP___ALGORITHM_SORT_H"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 1, "covered" : 0,};
var merged_data = [];

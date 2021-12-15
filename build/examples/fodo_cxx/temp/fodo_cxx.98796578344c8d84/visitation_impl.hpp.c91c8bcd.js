var data = {lines:[
{"lineNum":"    1","line":"//-----------------------------------------------------------------------------"},
{"lineNum":"    2","line":"// boost variant/detail/visitation_impl.hpp header file"},
{"lineNum":"    3","line":"// See http://www.boost.org for updates, documentation, and revision history."},
{"lineNum":"    4","line":"//-----------------------------------------------------------------------------"},
{"lineNum":"    5","line":"//"},
{"lineNum":"    6","line":"// Copyright (c) 2003"},
{"lineNum":"    7","line":"// Eric Friedman"},
{"lineNum":"    8","line":"//"},
{"lineNum":"    9","line":"// Distributed under the Boost Software License, Version 1.0. (See"},
{"lineNum":"   10","line":"// accompanying file LICENSE_1_0.txt or copy at"},
{"lineNum":"   11","line":"// http://www.boost.org/LICENSE_1_0.txt)"},
{"lineNum":"   12","line":""},
{"lineNum":"   13","line":"#ifndef BOOST_VARIANT_DETAIL_VISITATION_IMPL_HPP"},
{"lineNum":"   14","line":"#define BOOST_VARIANT_DETAIL_VISITATION_IMPL_HPP"},
{"lineNum":"   15","line":""},
{"lineNum":"   16","line":"#include <boost/config.hpp>"},
{"lineNum":"   17","line":""},
{"lineNum":"   18","line":"#include <boost/variant/detail/backup_holder.hpp>"},
{"lineNum":"   19","line":"#include <boost/variant/detail/cast_storage.hpp>"},
{"lineNum":"   20","line":"#include <boost/variant/detail/forced_return.hpp>"},
{"lineNum":"   21","line":"#include <boost/variant/variant_fwd.hpp> // for BOOST_VARIANT_DO_NOT_USE_VARIADIC_TEMPLATES"},
{"lineNum":"   22","line":""},
{"lineNum":"   23","line":"#include <boost/mpl/eval_if.hpp>"},
{"lineNum":"   24","line":"#include <boost/mpl/bool.hpp>"},
{"lineNum":"   25","line":"#include <boost/mpl/identity.hpp>"},
{"lineNum":"   26","line":"#include <boost/mpl/int.hpp>"},
{"lineNum":"   27","line":"#include <boost/mpl/next.hpp>"},
{"lineNum":"   28","line":"#include <boost/mpl/deref.hpp>"},
{"lineNum":"   29","line":"#include <boost/mpl/or.hpp>"},
{"lineNum":"   30","line":"#include <boost/preprocessor/cat.hpp>"},
{"lineNum":"   31","line":"#include <boost/preprocessor/inc.hpp>"},
{"lineNum":"   32","line":"#include <boost/preprocessor/repeat.hpp>"},
{"lineNum":"   33","line":"#include <boost/type_traits/is_same.hpp>"},
{"lineNum":"   34","line":"#include <boost/type_traits/has_nothrow_copy.hpp>"},
{"lineNum":"   35","line":"#include <boost/type_traits/is_nothrow_move_constructible.hpp>"},
{"lineNum":"   36","line":""},
{"lineNum":"   37","line":"#if BOOST_WORKAROUND(BOOST_MSVC, >= 1400)"},
{"lineNum":"   38","line":"# pragma warning (push)"},
{"lineNum":"   39","line":"# pragma warning (disable : 4702) //unreachable code"},
{"lineNum":"   40","line":"#endif"},
{"lineNum":"   41","line":""},
{"lineNum":"   42","line":"///////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   43","line":"// BOOST_VARIANT_VISITATION_UNROLLING_LIMIT"},
{"lineNum":"   44","line":"//"},
{"lineNum":"   45","line":"// Unrolls variant\'s visitation mechanism to reduce template instantiation"},
{"lineNum":"   46","line":"// and potentially increase runtime performance. (TODO: Investigate further.)"},
{"lineNum":"   47","line":"//"},
{"lineNum":"   48","line":"#if !defined(BOOST_VARIANT_VISITATION_UNROLLING_LIMIT)"},
{"lineNum":"   49","line":""},
{"lineNum":"   50","line":"#ifndef BOOST_VARIANT_DO_NOT_USE_VARIADIC_TEMPLATES"},
{"lineNum":"   51","line":"#   include <boost/mpl/limits/list.hpp>"},
{"lineNum":"   52","line":"#   define BOOST_VARIANT_VISITATION_UNROLLING_LIMIT   \\"},
{"lineNum":"   53","line":"        BOOST_MPL_LIMIT_LIST_SIZE"},
{"lineNum":"   54","line":"#else"},
{"lineNum":"   55","line":"#   define BOOST_VARIANT_VISITATION_UNROLLING_LIMIT   \\"},
{"lineNum":"   56","line":"        BOOST_VARIANT_LIMIT_TYPES"},
{"lineNum":"   57","line":"#endif"},
{"lineNum":"   58","line":""},
{"lineNum":"   59","line":"#endif"},
{"lineNum":"   60","line":""},
{"lineNum":"   61","line":"namespace boost {"},
{"lineNum":"   62","line":"namespace detail { namespace variant {"},
{"lineNum":"   63","line":""},
{"lineNum":"   64","line":"///////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   65","line":"// (detail) class apply_visitor_unrolled"},
{"lineNum":"   66","line":"//"},
{"lineNum":"   67","line":"// Tag type indicates when visitation_impl is unrolled."},
{"lineNum":"   68","line":"//"},
{"lineNum":"   69","line":"struct apply_visitor_unrolled {};"},
{"lineNum":"   70","line":""},
{"lineNum":"   71","line":"///////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   72","line":"// (detail) class template visitation_impl_step"},
{"lineNum":"   73","line":"//"},
{"lineNum":"   74","line":"// \"Never ending\" iterator range facilitates visitation_impl unrolling."},
{"lineNum":"   75","line":"//"},
{"lineNum":"   76","line":""},
{"lineNum":"   77","line":""},
{"lineNum":"   78","line":"template <typename Iter, typename LastIter>"},
{"lineNum":"   79","line":"struct visitation_impl_step"},
{"lineNum":"   80","line":"{"},
{"lineNum":"   81","line":"    typedef typename mpl::deref<Iter>::type type;"},
{"lineNum":"   82","line":""},
{"lineNum":"   83","line":"    typedef typename mpl::next<Iter>::type next_iter;"},
{"lineNum":"   84","line":"    typedef visitation_impl_step<"},
{"lineNum":"   85","line":"          next_iter, LastIter"},
{"lineNum":"   86","line":"        > next;"},
{"lineNum":"   87","line":"};"},
{"lineNum":"   88","line":""},
{"lineNum":"   89","line":"template <typename LastIter>"},
{"lineNum":"   90","line":"struct visitation_impl_step< LastIter,LastIter >"},
{"lineNum":"   91","line":"{"},
{"lineNum":"   92","line":"    typedef apply_visitor_unrolled type;"},
{"lineNum":"   93","line":"    typedef visitation_impl_step next;"},
{"lineNum":"   94","line":"};"},
{"lineNum":"   95","line":""},
{"lineNum":"   96","line":""},
{"lineNum":"   97","line":"///////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   98","line":"// (detail) function template visitation_impl_invoke"},
{"lineNum":"   99","line":"//"},
{"lineNum":"  100","line":"// Invokes the given visitor on the specified type in the given storage."},
{"lineNum":"  101","line":"//"},
{"lineNum":"  102","line":""},
{"lineNum":"  103","line":"template <typename Visitor, typename VoidPtrCV, typename T>"},
{"lineNum":"  104","line":"inline typename Visitor::result_type"},
{"lineNum":"  105","line":"visitation_impl_invoke_impl("},
{"lineNum":"  106","line":"      int, Visitor& visitor, VoidPtrCV storage, T*"},
{"lineNum":"  107","line":"    , mpl::true_// never_uses_backup"},
{"lineNum":"  108","line":"    )"},
{"lineNum":"  109","line":"{"},
{"lineNum":"  110","line":"    return visitor.internal_visit("},
{"lineNum":"  111","line":"          cast_storage<T>(storage), 1L"},
{"lineNum":"  112","line":"        );"},
{"lineNum":"  113","line":"}"},
{"lineNum":"  114","line":""},
{"lineNum":"  115","line":"template <typename Visitor, typename VoidPtrCV, typename T>"},
{"lineNum":"  116","line":"inline typename Visitor::result_type"},
{"lineNum":"  117","line":"visitation_impl_invoke_impl("},
{"lineNum":"  118","line":"      int internal_which, Visitor& visitor, VoidPtrCV storage, T*"},
{"lineNum":"  119","line":"    , mpl::false_// never_uses_backup"},
{"lineNum":"  120","line":"    )"},
{"lineNum":"  121","line":"{"},
{"lineNum":"  122","line":"    if (internal_which >= 0)"},
{"lineNum":"  123","line":"    {"},
{"lineNum":"  124","line":"        return visitor.internal_visit("},
{"lineNum":"  125","line":"              cast_storage<T>(storage), 1L"},
{"lineNum":"  126","line":"            );"},
{"lineNum":"  127","line":"    }"},
{"lineNum":"  128","line":"    else"},
{"lineNum":"  129","line":"    {"},
{"lineNum":"  130","line":"        return visitor.internal_visit("},
{"lineNum":"  131","line":"              cast_storage< backup_holder<T> >(storage), 1L"},
{"lineNum":"  132","line":"            );"},
{"lineNum":"  133","line":"    }"},
{"lineNum":"  134","line":"}"},
{"lineNum":"  135","line":""},
{"lineNum":"  136","line":"template <typename Visitor, typename VoidPtrCV, typename T, typename NoBackupFlag>"},
{"lineNum":"  137","line":"inline typename Visitor::result_type"},
{"lineNum":"  138","line":"visitation_impl_invoke("},
{"lineNum":"  139","line":"      int internal_which, Visitor& visitor, VoidPtrCV storage, T* t"},
{"lineNum":"  140","line":"    , NoBackupFlag"},
{"lineNum":"  141","line":"    , int"},
{"lineNum":"  142","line":"    )"},
{"lineNum":"  143","line":"{"},
{"lineNum":"  144","line":"    typedef typename mpl::or_<"},
{"lineNum":"  145","line":"          NoBackupFlag"},
{"lineNum":"  146","line":"        , is_nothrow_move_constructible<T>"},
{"lineNum":"  147","line":"        , has_nothrow_copy<T>"},
{"lineNum":"  148","line":"        >::type never_uses_backup;"},
{"lineNum":"  149","line":""},
{"lineNum":"  150","line":"    return (visitation_impl_invoke_impl)("},
{"lineNum":"  151","line":"          internal_which, visitor, storage, t"},
{"lineNum":"  152","line":"        , never_uses_backup()"},
{"lineNum":"  153","line":"        );"},
{"lineNum":"  154","line":"}"},
{"lineNum":"  155","line":""},
{"lineNum":"  156","line":"template <typename Visitor, typename VoidPtrCV, typename NBF>"},
{"lineNum":"  157","line":"inline typename Visitor::result_type"},
{"lineNum":"  158","line":"visitation_impl_invoke(int, Visitor&, VoidPtrCV, apply_visitor_unrolled*, NBF, long)"},
{"lineNum":"  159","line":"{"},
{"lineNum":"  160","line":"    // should never be here at runtime!"},
{"lineNum":"  161","line":"    typedef typename Visitor::result_type result_type;"},
{"lineNum":"  162","line":"    return ::boost::detail::variant::forced_return< result_type >();"},
{"lineNum":"  163","line":"}"},
{"lineNum":"  164","line":""},
{"lineNum":"  165","line":"///////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"  166","line":"// (detail) function template visitation_impl"},
{"lineNum":"  167","line":"//"},
{"lineNum":"  168","line":"// Invokes the given visitor on the type in the given variant storage."},
{"lineNum":"  169","line":"//"},
{"lineNum":"  170","line":""},
{"lineNum":"  171","line":"template <"},
{"lineNum":"  172","line":"      typename W, typename S"},
{"lineNum":"  173","line":"    , typename Visitor, typename VPCV"},
{"lineNum":"  174","line":"    , typename NBF"},
{"lineNum":"  175","line":"    >"},
{"lineNum":"  176","line":"inline typename Visitor::result_type"},
{"lineNum":"  177","line":"visitation_impl("},
{"lineNum":"  178","line":"      int, int, Visitor&, VPCV"},
{"lineNum":"  179","line":"    , mpl::true_ // is_apply_visitor_unrolled"},
{"lineNum":"  180","line":"    , NBF, W* = 0, S* = 0"},
{"lineNum":"  181","line":"    )"},
{"lineNum":"  182","line":"{"},
{"lineNum":"  183","line":"    // should never be here at runtime!"},
{"lineNum":"  184","line":"    typedef typename Visitor::result_type result_type;"},
{"lineNum":"  185","line":"    return ::boost::detail::variant::forced_return< result_type >();"},
{"lineNum":"  186","line":"}"},
{"lineNum":"  187","line":""},
{"lineNum":"  188","line":"template <"},
{"lineNum":"  189","line":"      typename Which, typename step0"},
{"lineNum":"  190","line":"    , typename Visitor, typename VoidPtrCV"},
{"lineNum":"  191","line":"    , typename NoBackupFlag"},
{"lineNum":"  192","line":"    >"},
{"lineNum":"  193","line":"inline typename Visitor::result_type"},
{"lineNum":"  194","line":"visitation_impl("},
{"lineNum":"  195","line":"      const int internal_which, const int logical_which"},
{"lineNum":"  196","line":"    , Visitor& visitor, VoidPtrCV storage"},
{"lineNum":"  197","line":"    , mpl::false_ // is_apply_visitor_unrolled"},
{"lineNum":"  198","line":"    , NoBackupFlag no_backup_flag"},
{"lineNum":"  199","line":"    , Which* = 0, step0* = 0"},
{"lineNum":"  200","line":"    )"},
{"lineNum":"  201","line":"{","class":"lineNoCov","hits":"0","possible_hits":"18",},
{"lineNum":"  202","line":"    // Typedef apply_visitor_unrolled steps and associated types..."},
{"lineNum":"  203","line":"#   define BOOST_VARIANT_AUX_APPLY_VISITOR_STEP_TYPEDEF(z, N, _) \\"},
{"lineNum":"  204","line":"    typedef typename BOOST_PP_CAT(step,N)::type BOOST_PP_CAT(T,N); \\"},
{"lineNum":"  205","line":"    typedef typename BOOST_PP_CAT(step,N)::next \\"},
{"lineNum":"  206","line":"        BOOST_PP_CAT(step, BOOST_PP_INC(N)); \\"},
{"lineNum":"  207","line":"    /**/"},
{"lineNum":"  208","line":""},
{"lineNum":"  209","line":"    BOOST_PP_REPEAT("},
{"lineNum":"  210","line":"          BOOST_VARIANT_VISITATION_UNROLLING_LIMIT"},
{"lineNum":"  211","line":"        , BOOST_VARIANT_AUX_APPLY_VISITOR_STEP_TYPEDEF"},
{"lineNum":"  212","line":"        , _"},
{"lineNum":"  213","line":"        )"},
{"lineNum":"  214","line":""},
{"lineNum":"  215","line":"#   undef BOOST_VARIANT_AUX_APPLY_VISITOR_STEP_TYPEDEF"},
{"lineNum":"  216","line":""},
{"lineNum":"  217","line":"    // ...switch on the target which-index value..."},
{"lineNum":"  218","line":"    switch (logical_which)","class":"lineNoCov","hits":"0","possible_hits":"18",},
{"lineNum":"  219","line":"    {"},
{"lineNum":"  220","line":""},
{"lineNum":"  221","line":"    // ...applying the appropriate case:"},
{"lineNum":"  222","line":"#   define BOOST_VARIANT_AUX_APPLY_VISITOR_STEP_CASE(z, N, _) \\"},
{"lineNum":"  223","line":"    case (Which::value + (N)): \\"},
{"lineNum":"  224","line":"        return (visitation_impl_invoke)( \\"},
{"lineNum":"  225","line":"              internal_which, visitor, storage \\"},
{"lineNum":"  226","line":"            , static_cast<BOOST_PP_CAT(T,N)*>(0) \\"},
{"lineNum":"  227","line":"            , no_backup_flag, 1L \\"},
{"lineNum":"  228","line":"            ); \\"},
{"lineNum":"  229","line":"    /**/"},
{"lineNum":"  230","line":""},
{"lineNum":"  231","line":"    BOOST_PP_REPEAT("},
{"lineNum":"  232","line":"          BOOST_VARIANT_VISITATION_UNROLLING_LIMIT"},
{"lineNum":"  233","line":"        , BOOST_VARIANT_AUX_APPLY_VISITOR_STEP_CASE"},
{"lineNum":"  234","line":"        , _"},
{"lineNum":"  235","line":"        )"},
{"lineNum":"  236","line":""},
{"lineNum":"  237","line":"#   undef BOOST_VARIANT_AUX_APPLY_VISITOR_STEP_CASE"},
{"lineNum":"  238","line":""},
{"lineNum":"  239","line":"    default: break;"},
{"lineNum":"  240","line":"    }"},
{"lineNum":"  241","line":""},
{"lineNum":"  242","line":"    // If not handled in this iteration, continue unrolling:"},
{"lineNum":"  243","line":"    typedef mpl::int_<"},
{"lineNum":"  244","line":"          Which::value + (BOOST_VARIANT_VISITATION_UNROLLING_LIMIT)"},
{"lineNum":"  245","line":"        > next_which;"},
{"lineNum":"  246","line":""},
{"lineNum":"  247","line":"    typedef BOOST_PP_CAT(step, BOOST_VARIANT_VISITATION_UNROLLING_LIMIT)"},
{"lineNum":"  248","line":"        next_step;"},
{"lineNum":"  249","line":""},
{"lineNum":"  250","line":"    typedef typename next_step::type next_type;"},
{"lineNum":"  251","line":"    typedef typename is_same< next_type,apply_visitor_unrolled >::type"},
{"lineNum":"  252","line":"        is_apply_visitor_unrolled;"},
{"lineNum":"  253","line":""},
{"lineNum":"  254","line":"    return detail::variant::visitation_impl("},
{"lineNum":"  255","line":"          internal_which, logical_which"},
{"lineNum":"  256","line":"        , visitor, storage"},
{"lineNum":"  257","line":"        , is_apply_visitor_unrolled()"},
{"lineNum":"  258","line":"        , no_backup_flag"},
{"lineNum":"  259","line":"        , static_cast<next_which*>(0), static_cast<next_step*>(0)"},
{"lineNum":"  260","line":"        );"},
{"lineNum":"  261","line":"}","class":"lineNoCov","hits":"0","possible_hits":"22",},
{"lineNum":"  262","line":""},
{"lineNum":"  263","line":"}} // namespace detail::variant"},
{"lineNum":"  264","line":"} // namespace boost"},
{"lineNum":"  265","line":""},
{"lineNum":"  266","line":"#if BOOST_WORKAROUND(BOOST_MSVC, >= 1400)"},
{"lineNum":"  267","line":"# pragma warning(pop)"},
{"lineNum":"  268","line":"#endif"},
{"lineNum":"  269","line":""},
{"lineNum":"  270","line":"#endif // BOOST_VARIANT_DETAIL_VISITATION_IMPL_HPP"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 3, "covered" : 0,};
var merged_data = [];

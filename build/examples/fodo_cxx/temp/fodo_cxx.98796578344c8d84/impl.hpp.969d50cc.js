var data = {lines:[
{"lineNum":"    1","line":"///////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"    2","line":"/// \\file impl.hpp"},
{"lineNum":"    3","line":"/// Contains definition of transform<> and transform_impl<> helpers."},
{"lineNum":"    4","line":"//"},
{"lineNum":"    5","line":"//  Copyright 2008 Eric Niebler. Distributed under the Boost"},
{"lineNum":"    6","line":"//  Software License, Version 1.0. (See accompanying file"},
{"lineNum":"    7","line":"//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)"},
{"lineNum":"    8","line":""},
{"lineNum":"    9","line":"#ifndef BOOST_PROTO_TRANSFORM_IMPL_HPP_EAN_04_03_2008"},
{"lineNum":"   10","line":"#define BOOST_PROTO_TRANSFORM_IMPL_HPP_EAN_04_03_2008"},
{"lineNum":"   11","line":""},
{"lineNum":"   12","line":"#include <boost/config.hpp>"},
{"lineNum":"   13","line":"#include <boost/mpl/bool.hpp>"},
{"lineNum":"   14","line":"#include <boost/type_traits/add_const.hpp>"},
{"lineNum":"   15","line":"#include <boost/type_traits/add_reference.hpp>"},
{"lineNum":"   16","line":"#include <boost/proto/proto_fwd.hpp>"},
{"lineNum":"   17","line":"#include <boost/proto/detail/any.hpp>"},
{"lineNum":"   18","line":"#include <boost/proto/detail/static_const.hpp>"},
{"lineNum":"   19","line":""},
{"lineNum":"   20","line":"#if defined(_MSC_VER)"},
{"lineNum":"   21","line":"# pragma warning(push)"},
{"lineNum":"   22","line":"# pragma warning(disable : 4714) // function \'xxx\' marked as __forceinline not inlined"},
{"lineNum":"   23","line":"#endif"},
{"lineNum":"   24","line":""},
{"lineNum":"   25","line":"namespace boost { namespace proto"},
{"lineNum":"   26","line":"{"},
{"lineNum":"   27","line":"    namespace envns_"},
{"lineNum":"   28","line":"    {"},
{"lineNum":"   29","line":"        ////////////////////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   30","line":"        struct key_not_found"},
{"lineNum":"   31","line":"        {};"},
{"lineNum":"   32","line":""},
{"lineNum":"   33","line":"        ////////////////////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   34","line":"        // empty_env"},
{"lineNum":"   35","line":"        struct empty_env"},
{"lineNum":"   36","line":"        {"},
{"lineNum":"   37","line":"            typedef void proto_environment_;"},
{"lineNum":"   38","line":""},
{"lineNum":"   39","line":"            template<typename OtherTag, typename OtherValue = key_not_found>"},
{"lineNum":"   40","line":"            struct lookup"},
{"lineNum":"   41","line":"            {"},
{"lineNum":"   42","line":"                typedef OtherValue type;"},
{"lineNum":"   43","line":"                typedef"},
{"lineNum":"   44","line":"                    typename add_reference<typename add_const<OtherValue>::type>::type"},
{"lineNum":"   45","line":"                const_reference;"},
{"lineNum":"   46","line":"            };"},
{"lineNum":"   47","line":""},
{"lineNum":"   48","line":"            key_not_found operator[](detail::any) const"},
{"lineNum":"   49","line":"            {"},
{"lineNum":"   50","line":"                return key_not_found();"},
{"lineNum":"   51","line":"            }"},
{"lineNum":"   52","line":""},
{"lineNum":"   53","line":"            template<typename T>"},
{"lineNum":"   54","line":"            T const &at(detail::any, T const &t) const"},
{"lineNum":"   55","line":"            {"},
{"lineNum":"   56","line":"                return t;"},
{"lineNum":"   57","line":"            }"},
{"lineNum":"   58","line":"        };"},
{"lineNum":"   59","line":"    }"},
{"lineNum":"   60","line":""},
{"lineNum":"   61","line":"    ////////////////////////////////////////////////////////////////////////////////////////////"},
{"lineNum":"   62","line":"    // is_env"},
{"lineNum":"   63","line":"    template<typename T, typename Void>"},
{"lineNum":"   64","line":"    struct is_env"},
{"lineNum":"   65","line":"      : mpl::false_"},
{"lineNum":"   66","line":"    {};"},
{"lineNum":"   67","line":""},
{"lineNum":"   68","line":"    template<typename T>"},
{"lineNum":"   69","line":"    struct is_env<T, typename T::proto_environment_>"},
{"lineNum":"   70","line":"      : mpl::true_"},
{"lineNum":"   71","line":"    {};"},
{"lineNum":"   72","line":""},
{"lineNum":"   73","line":"    template<typename T>"},
{"lineNum":"   74","line":"    struct is_env<T &, void>"},
{"lineNum":"   75","line":"      : is_env<T>"},
{"lineNum":"   76","line":"    {};"},
{"lineNum":"   77","line":""},
{"lineNum":"   78","line":"#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES"},
{"lineNum":"   79","line":""},
{"lineNum":"   80","line":"    /// INTERNAL ONLY"},
{"lineNum":"   81","line":"    ///"},
{"lineNum":"   82","line":"    #define BOOST_PROTO_TRANSFORM_(PrimitiveTransform, X)                                                       \\"},
{"lineNum":"   83","line":"    BOOST_PROTO_CALLABLE()                                                                                      \\"},
{"lineNum":"   84","line":"    typedef X proto_is_transform_;                                                                              \\"},
{"lineNum":"   85","line":"    typedef PrimitiveTransform transform_type;                                                                  \\"},
{"lineNum":"   86","line":"                                                                                                                \\"},
{"lineNum":"   87","line":"    template<typename Sig>                                                                                      \\"},
{"lineNum":"   88","line":"    struct result                                                                                               \\"},
{"lineNum":"   89","line":"    {                                                                                                           \\"},
{"lineNum":"   90","line":"        typedef typename boost::proto::detail::apply_transform<Sig>::result_type type;                          \\"},
{"lineNum":"   91","line":"    };                                                                                                          \\"},
{"lineNum":"   92","line":"                                                                                                                \\"},
{"lineNum":"   93","line":"    template<typename Expr>                                                                                     \\"},
{"lineNum":"   94","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"   95","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr &)>::result_type                         \\"},
{"lineNum":"   96","line":"    operator ()(Expr &e) const                                                                                  \\"},
{"lineNum":"   97","line":"    {                                                                                                           \\"},
{"lineNum":"   98","line":"        boost::proto::empty_state s = 0;                                                                        \\"},
{"lineNum":"   99","line":"        boost::proto::empty_env d;                                                                              \\"},
{"lineNum":"  100","line":"        return boost::proto::detail::apply_transform<transform_type(Expr &)>()(e, s, d);                        \\"},
{"lineNum":"  101","line":"    }                                                                                                           \\"},
{"lineNum":"  102","line":"                                                                                                                \\"},
{"lineNum":"  103","line":"    template<typename Expr>                                                                                     \\"},
{"lineNum":"  104","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"  105","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr const &)>::result_type                   \\"},
{"lineNum":"  106","line":"    operator ()(Expr const &e) const                                                                            \\"},
{"lineNum":"  107","line":"    {                                                                                                           \\"},
{"lineNum":"  108","line":"        boost::proto::empty_state s = 0;                                                                        \\"},
{"lineNum":"  109","line":"        boost::proto::empty_env d;                                                                              \\"},
{"lineNum":"  110","line":"        return boost::proto::detail::apply_transform<transform_type(Expr const &)>()(e, s, d);                  \\"},
{"lineNum":"  111","line":"    }                                                                                                           \\"},
{"lineNum":"  112","line":"                                                                                                                \\"},
{"lineNum":"  113","line":"    template<typename Expr, typename State>                                                                     \\"},
{"lineNum":"  114","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"  115","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr &, State &)>::result_type                \\"},
{"lineNum":"  116","line":"    operator ()(Expr &e, State &s) const                                                                        \\"},
{"lineNum":"  117","line":"    {                                                                                                           \\"},
{"lineNum":"  118","line":"        boost::proto::empty_env d;                                                                              \\"},
{"lineNum":"  119","line":"        return boost::proto::detail::apply_transform<transform_type(Expr &, State &)>()(e, s, d);               \\"},
{"lineNum":"  120","line":"    }                                                                                                           \\"},
{"lineNum":"  121","line":"                                                                                                                \\"},
{"lineNum":"  122","line":"    template<typename Expr, typename State>                                                                     \\"},
{"lineNum":"  123","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"  124","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr const &, State &)>::result_type          \\"},
{"lineNum":"  125","line":"    operator ()(Expr const &e, State &s) const                                                                  \\"},
{"lineNum":"  126","line":"    {                                                                                                           \\"},
{"lineNum":"  127","line":"        boost::proto::empty_env d;                                                                              \\"},
{"lineNum":"  128","line":"        return boost::proto::detail::apply_transform<transform_type(Expr const &, State &)>()(e, s, d);         \\"},
{"lineNum":"  129","line":"    }                                                                                                           \\"},
{"lineNum":"  130","line":"                                                                                                                \\"},
{"lineNum":"  131","line":"    template<typename Expr, typename State>                                                                     \\"},
{"lineNum":"  132","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"  133","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr &, State const &)>::result_type          \\"},
{"lineNum":"  134","line":"    operator ()(Expr &e, State const &s) const                                                                  \\"},
{"lineNum":"  135","line":"    {                                                                                                           \\"},
{"lineNum":"  136","line":"        boost::proto::empty_env d;                                                                              \\"},
{"lineNum":"  137","line":"        return boost::proto::detail::apply_transform<transform_type(Expr &, State const &)>()(e, s, d);         \\"},
{"lineNum":"  138","line":"    }                                                                                                           \\"},
{"lineNum":"  139","line":"                                                                                                                \\"},
{"lineNum":"  140","line":"    template<typename Expr, typename State>                                                                     \\"},
{"lineNum":"  141","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"  142","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr const &, State const &)>::result_type    \\"},
{"lineNum":"  143","line":"    operator ()(Expr const &e, State const &s) const                                                            \\"},
{"lineNum":"  144","line":"    {                                                                                                           \\"},
{"lineNum":"  145","line":"        boost::proto::empty_env d;                                                                              \\"},
{"lineNum":"  146","line":"        return boost::proto::detail::apply_transform<transform_type(Expr const &, State const &)>()(e, s, d);   \\"},
{"lineNum":"  147","line":"    }                                                                                                           \\"},
{"lineNum":"  148","line":"                                                                                                                \\"},
{"lineNum":"  149","line":"    template<typename Expr, typename State, typename Data>                                                      \\"},
{"lineNum":"  150","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"  151","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr &, State &, Data &)>::result_type        \\"},
{"lineNum":"  152","line":"    operator ()(Expr &e, State &s, Data &d) const                                                               \\"},
{"lineNum":"  153","line":"    {                                                                                                           \\"},
{"lineNum":"  154","line":"        return boost::proto::detail::apply_transform<transform_type(Expr &, State &, Data &)>()(e, s, d);       \\"},
{"lineNum":"  155","line":"    }                                                                                                           \\"},
{"lineNum":"  156","line":"                                                                                                                \\"},
{"lineNum":"  157","line":"    template<typename Expr, typename State, typename Data>                                                      \\"},
{"lineNum":"  158","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"  159","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr const &, State &, Data &)>::result_type  \\"},
{"lineNum":"  160","line":"    operator ()(Expr const &e, State &s, Data &d) const                                                         \\"},
{"lineNum":"  161","line":"    {                                                                                                           \\"},
{"lineNum":"  162","line":"        return boost::proto::detail::apply_transform<transform_type(Expr const &, State &, Data &)>()(e, s, d); \\"},
{"lineNum":"  163","line":"    }                                                                                                           \\"},
{"lineNum":"  164","line":"                                                                                                                \\"},
{"lineNum":"  165","line":"    template<typename Expr, typename State, typename Data>                                                      \\"},
{"lineNum":"  166","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"  167","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr &, State const &, Data &)>::result_type  \\"},
{"lineNum":"  168","line":"    operator ()(Expr &e, State const &s, Data &d) const                                                         \\"},
{"lineNum":"  169","line":"    {                                                                                                           \\"},
{"lineNum":"  170","line":"        return boost::proto::detail::apply_transform<transform_type(Expr &, State const &, Data &)>()(e, s, d); \\"},
{"lineNum":"  171","line":"    }                                                                                                           \\"},
{"lineNum":"  172","line":"                                                                                                                \\"},
{"lineNum":"  173","line":"    template<typename Expr, typename State, typename Data>                                                      \\"},
{"lineNum":"  174","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"  175","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr const &, State const &, Data &)>::result_type  \\"},
{"lineNum":"  176","line":"    operator ()(Expr const &e, State const &s, Data &d) const                                                   \\"},
{"lineNum":"  177","line":"    {                                                                                                           \\"},
{"lineNum":"  178","line":"        return boost::proto::detail::apply_transform<transform_type(Expr const &, State const &, Data &)>()(e, s, d); \\"},
{"lineNum":"  179","line":"    }                                                                                                           \\"},
{"lineNum":"  180","line":"    /**/"},
{"lineNum":"  181","line":""},
{"lineNum":"  182","line":"#else"},
{"lineNum":"  183","line":""},
{"lineNum":"  184","line":"    /// INTERNAL ONLY"},
{"lineNum":"  185","line":"    ///"},
{"lineNum":"  186","line":"    #define BOOST_PROTO_TRANSFORM_(PrimitiveTransform, X)                                                       \\"},
{"lineNum":"  187","line":"    BOOST_PROTO_CALLABLE()                                                                                      \\"},
{"lineNum":"  188","line":"    typedef X proto_is_transform_;                                                                              \\"},
{"lineNum":"  189","line":"    typedef PrimitiveTransform transform_type;                                                                  \\"},
{"lineNum":"  190","line":"                                                                                                                \\"},
{"lineNum":"  191","line":"    template<typename Sig>                                                                                      \\"},
{"lineNum":"  192","line":"    struct result                                                                                               \\"},
{"lineNum":"  193","line":"    {                                                                                                           \\"},
{"lineNum":"  194","line":"        typedef typename boost::proto::detail::apply_transform<Sig>::result_type type;                          \\"},
{"lineNum":"  195","line":"    };                                                                                                          \\"},
{"lineNum":"  196","line":"                                                                                                                \\"},
{"lineNum":"  197","line":"    template<typename Expr>                                                                                     \\"},
{"lineNum":"  198","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"  199","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr const &)>::result_type                   \\"},
{"lineNum":"  200","line":"    operator ()(Expr &&e) const                                                                                 \\"},
{"lineNum":"  201","line":"    {                                                                                                           \\"},
{"lineNum":"  202","line":"        boost::proto::empty_state s = 0;                                                                        \\"},
{"lineNum":"  203","line":"        boost::proto::empty_env d;                                                                              \\"},
{"lineNum":"  204","line":"        return boost::proto::detail::apply_transform<transform_type(Expr const &)>()(e, s, d);                  \\"},
{"lineNum":"  205","line":"    }                                                                                                           \\"},
{"lineNum":"  206","line":"                                                                                                                \\"},
{"lineNum":"  207","line":"    template<typename Expr, typename State>                                                                     \\"},
{"lineNum":"  208","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"  209","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr const &, State const &)>::result_type    \\"},
{"lineNum":"  210","line":"    operator ()(Expr &&e, State &&s) const                                                                      \\"},
{"lineNum":"  211","line":"    {                                                                                                           \\"},
{"lineNum":"  212","line":"        boost::proto::empty_env d;                                                                              \\"},
{"lineNum":"  213","line":"        return boost::proto::detail::apply_transform<transform_type(Expr const &, State const &)>()(e, s, d);   \\"},
{"lineNum":"  214","line":"    }                                                                                                           \\"},
{"lineNum":"  215","line":"                                                                                                                \\"},
{"lineNum":"  216","line":"    template<typename Expr, typename State, typename Data>                                                      \\"},
{"lineNum":"  217","line":"    BOOST_FORCEINLINE                                                                                           \\"},
{"lineNum":"  218","line":"    typename boost::proto::detail::apply_transform<transform_type(Expr const &, State const &, Data const &)>::result_type \\"},
{"lineNum":"  219","line":"    operator ()(Expr &&e, State &&s, Data &&d) const                                                            \\"},
{"lineNum":"  220","line":"    {                                                                                                           \\"},
{"lineNum":"  221","line":"        return boost::proto::detail::apply_transform<transform_type(Expr const &, State const &, Data const &)>()(e, s, d); \\"},
{"lineNum":"  222","line":"    }                                                                                                           \\"},
{"lineNum":"  223","line":"    /**/"},
{"lineNum":"  224","line":""},
{"lineNum":"  225","line":"#endif"},
{"lineNum":"  226","line":""},
{"lineNum":"  227","line":"    #define BOOST_PROTO_TRANSFORM(PrimitiveTransform)                                                           \\"},
{"lineNum":"  228","line":"        BOOST_PROTO_TRANSFORM_(PrimitiveTransform, void)                                                        \\"},
{"lineNum":"  229","line":"        /**/"},
{"lineNum":"  230","line":""},
{"lineNum":"  231","line":"    namespace detail"},
{"lineNum":"  232","line":"    {"},
{"lineNum":"  233","line":"        template<typename Sig>"},
{"lineNum":"  234","line":"        struct apply_transform;"},
{"lineNum":"  235","line":""},
{"lineNum":"  236","line":"        template<typename PrimitiveTransform, typename Expr>"},
{"lineNum":"  237","line":"        struct apply_transform<PrimitiveTransform(Expr)>"},
{"lineNum":"  238","line":"          : PrimitiveTransform::template impl<Expr, empty_state, empty_env>"},
{"lineNum":"  239","line":"        {};"},
{"lineNum":"  240","line":""},
{"lineNum":"  241","line":"        template<typename PrimitiveTransform, typename Expr, typename State>"},
{"lineNum":"  242","line":"        struct apply_transform<PrimitiveTransform(Expr, State)>"},
{"lineNum":"  243","line":"          : PrimitiveTransform::template impl<Expr, State, empty_env>"},
{"lineNum":"  244","line":"        {};"},
{"lineNum":"  245","line":""},
{"lineNum":"  246","line":"        template<typename PrimitiveTransform, typename Expr, typename State, typename Data>"},
{"lineNum":"  247","line":"        struct apply_transform<PrimitiveTransform(Expr, State, Data)>"},
{"lineNum":"  248","line":"          : PrimitiveTransform::template impl<Expr, State, Data>"},
{"lineNum":"  249","line":"        {};"},
{"lineNum":"  250","line":"    }"},
{"lineNum":"  251","line":""},
{"lineNum":"  252","line":"    template<typename PrimitiveTransform, typename X>"},
{"lineNum":"  253","line":"    struct transform"},
{"lineNum":"  254","line":"    {"},
{"lineNum":"  255","line":"        BOOST_PROTO_TRANSFORM_(PrimitiveTransform, X)","class":"lineNoCov","hits":"0","possible_hits":"29",},
{"lineNum":"  256","line":"    };"},
{"lineNum":"  257","line":""},
{"lineNum":"  258","line":"    template<typename Expr, typename State, typename Data>"},
{"lineNum":"  259","line":"    struct transform_impl"},
{"lineNum":"  260","line":"    {"},
{"lineNum":"  261","line":"        typedef Expr const expr;"},
{"lineNum":"  262","line":"        typedef Expr const &expr_param;"},
{"lineNum":"  263","line":"        typedef State const state;"},
{"lineNum":"  264","line":"        typedef State const &state_param;"},
{"lineNum":"  265","line":"        typedef Data const data;"},
{"lineNum":"  266","line":"        typedef Data const &data_param;"},
{"lineNum":"  267","line":"    };"},
{"lineNum":"  268","line":""},
{"lineNum":"  269","line":"    template<typename Expr, typename State, typename Data>"},
{"lineNum":"  270","line":"    struct transform_impl<Expr &, State, Data>"},
{"lineNum":"  271","line":"    {"},
{"lineNum":"  272","line":"        typedef Expr expr;"},
{"lineNum":"  273","line":"        typedef Expr &expr_param;"},
{"lineNum":"  274","line":"        typedef State const state;"},
{"lineNum":"  275","line":"        typedef State const &state_param;"},
{"lineNum":"  276","line":"        typedef Data const data;"},
{"lineNum":"  277","line":"        typedef Data const &data_param;"},
{"lineNum":"  278","line":"    };"},
{"lineNum":"  279","line":""},
{"lineNum":"  280","line":"    template<typename Expr, typename State, typename Data>"},
{"lineNum":"  281","line":"    struct transform_impl<Expr, State &, Data>"},
{"lineNum":"  282","line":"    {"},
{"lineNum":"  283","line":"        typedef Expr const expr;"},
{"lineNum":"  284","line":"        typedef Expr const &expr_param;"},
{"lineNum":"  285","line":"        typedef State state;"},
{"lineNum":"  286","line":"        typedef State &state_param;"},
{"lineNum":"  287","line":"        typedef Data const data;"},
{"lineNum":"  288","line":"        typedef Data const &data_param;"},
{"lineNum":"  289","line":"    };"},
{"lineNum":"  290","line":""},
{"lineNum":"  291","line":"    template<typename Expr, typename State, typename Data>"},
{"lineNum":"  292","line":"    struct transform_impl<Expr, State, Data &>"},
{"lineNum":"  293","line":"    {"},
{"lineNum":"  294","line":"        typedef Expr const expr;"},
{"lineNum":"  295","line":"        typedef Expr const &expr_param;"},
{"lineNum":"  296","line":"        typedef State const state;"},
{"lineNum":"  297","line":"        typedef State const &state_param;"},
{"lineNum":"  298","line":"        typedef Data data;"},
{"lineNum":"  299","line":"        typedef Data &data_param;"},
{"lineNum":"  300","line":"    };"},
{"lineNum":"  301","line":""},
{"lineNum":"  302","line":"    template<typename Expr, typename State, typename Data>"},
{"lineNum":"  303","line":"    struct transform_impl<Expr &, State &, Data>"},
{"lineNum":"  304","line":"    {"},
{"lineNum":"  305","line":"        typedef Expr expr;"},
{"lineNum":"  306","line":"        typedef Expr &expr_param;"},
{"lineNum":"  307","line":"        typedef State state;"},
{"lineNum":"  308","line":"        typedef State &state_param;"},
{"lineNum":"  309","line":"        typedef Data const data;"},
{"lineNum":"  310","line":"        typedef Data const &data_param;"},
{"lineNum":"  311","line":"    };"},
{"lineNum":"  312","line":""},
{"lineNum":"  313","line":"    template<typename Expr, typename State, typename Data>"},
{"lineNum":"  314","line":"    struct transform_impl<Expr &, State, Data &>"},
{"lineNum":"  315","line":"    {"},
{"lineNum":"  316","line":"        typedef Expr expr;"},
{"lineNum":"  317","line":"        typedef Expr &expr_param;"},
{"lineNum":"  318","line":"        typedef State const state;"},
{"lineNum":"  319","line":"        typedef State const &state_param;"},
{"lineNum":"  320","line":"        typedef Data data;"},
{"lineNum":"  321","line":"        typedef Data &data_param;"},
{"lineNum":"  322","line":"    };"},
{"lineNum":"  323","line":""},
{"lineNum":"  324","line":"    template<typename Expr, typename State, typename Data>"},
{"lineNum":"  325","line":"    struct transform_impl<Expr, State &, Data &>"},
{"lineNum":"  326","line":"    {"},
{"lineNum":"  327","line":"        typedef Expr const expr;"},
{"lineNum":"  328","line":"        typedef Expr const &expr_param;"},
{"lineNum":"  329","line":"        typedef State state;"},
{"lineNum":"  330","line":"        typedef State &state_param;"},
{"lineNum":"  331","line":"        typedef Data data;"},
{"lineNum":"  332","line":"        typedef Data &data_param;"},
{"lineNum":"  333","line":"    };"},
{"lineNum":"  334","line":""},
{"lineNum":"  335","line":"    template<typename Expr, typename State, typename Data>"},
{"lineNum":"  336","line":"    struct transform_impl<Expr &, State &, Data &>"},
{"lineNum":"  337","line":"    {"},
{"lineNum":"  338","line":"        typedef Expr expr;"},
{"lineNum":"  339","line":"        typedef Expr &expr_param;"},
{"lineNum":"  340","line":"        typedef State state;"},
{"lineNum":"  341","line":"        typedef State &state_param;"},
{"lineNum":"  342","line":"        typedef Data data;"},
{"lineNum":"  343","line":"        typedef Data &data_param;"},
{"lineNum":"  344","line":"    };"},
{"lineNum":"  345","line":""},
{"lineNum":"  346","line":"}} // namespace boost::proto"},
{"lineNum":"  347","line":""},
{"lineNum":"  348","line":"#if defined(_MSC_VER)"},
{"lineNum":"  349","line":"# pragma warning(pop)"},
{"lineNum":"  350","line":"#endif"},
{"lineNum":"  351","line":""},
{"lineNum":"  352","line":"#endif"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 1, "covered" : 0,};
var merged_data = [];

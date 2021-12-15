var data = {lines:[
{"lineNum":"    1","line":"// -*- C++ -*-"},
{"lineNum":"    2","line":"//===--------------------------- thread -----------------------------------===//"},
{"lineNum":"    3","line":"//"},
{"lineNum":"    4","line":"// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions."},
{"lineNum":"    5","line":"// See https://llvm.org/LICENSE.txt for license information."},
{"lineNum":"    6","line":"// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception"},
{"lineNum":"    7","line":"//"},
{"lineNum":"    8","line":"//===----------------------------------------------------------------------===//"},
{"lineNum":"    9","line":""},
{"lineNum":"   10","line":"#ifndef _LIBCPP_THREAD"},
{"lineNum":"   11","line":"#define _LIBCPP_THREAD"},
{"lineNum":"   12","line":""},
{"lineNum":"   13","line":"/*"},
{"lineNum":"   14","line":""},
{"lineNum":"   15","line":"    thread synopsis"},
{"lineNum":"   16","line":""},
{"lineNum":"   17","line":"namespace std"},
{"lineNum":"   18","line":"{"},
{"lineNum":"   19","line":""},
{"lineNum":"   20","line":"class thread"},
{"lineNum":"   21","line":"{"},
{"lineNum":"   22","line":"public:"},
{"lineNum":"   23","line":"    class id;"},
{"lineNum":"   24","line":"    typedef pthread_t native_handle_type;"},
{"lineNum":"   25","line":""},
{"lineNum":"   26","line":"    thread() noexcept;"},
{"lineNum":"   27","line":"    template <class F, class ...Args> explicit thread(F&& f, Args&&... args);"},
{"lineNum":"   28","line":"    ~thread();"},
{"lineNum":"   29","line":""},
{"lineNum":"   30","line":"    thread(const thread&) = delete;"},
{"lineNum":"   31","line":"    thread(thread&& t) noexcept;"},
{"lineNum":"   32","line":""},
{"lineNum":"   33","line":"    thread& operator=(const thread&) = delete;"},
{"lineNum":"   34","line":"    thread& operator=(thread&& t) noexcept;"},
{"lineNum":"   35","line":""},
{"lineNum":"   36","line":"    void swap(thread& t) noexcept;"},
{"lineNum":"   37","line":""},
{"lineNum":"   38","line":"    bool joinable() const noexcept;"},
{"lineNum":"   39","line":"    void join();"},
{"lineNum":"   40","line":"    void detach();"},
{"lineNum":"   41","line":"    id get_id() const noexcept;"},
{"lineNum":"   42","line":"    native_handle_type native_handle();"},
{"lineNum":"   43","line":""},
{"lineNum":"   44","line":"    static unsigned hardware_concurrency() noexcept;"},
{"lineNum":"   45","line":"};"},
{"lineNum":"   46","line":""},
{"lineNum":"   47","line":"void swap(thread& x, thread& y) noexcept;"},
{"lineNum":"   48","line":""},
{"lineNum":"   49","line":"class thread::id"},
{"lineNum":"   50","line":"{"},
{"lineNum":"   51","line":"public:"},
{"lineNum":"   52","line":"    id() noexcept;"},
{"lineNum":"   53","line":"};"},
{"lineNum":"   54","line":""},
{"lineNum":"   55","line":"bool operator==(thread::id x, thread::id y) noexcept;"},
{"lineNum":"   56","line":"bool operator!=(thread::id x, thread::id y) noexcept;"},
{"lineNum":"   57","line":"bool operator< (thread::id x, thread::id y) noexcept;"},
{"lineNum":"   58","line":"bool operator<=(thread::id x, thread::id y) noexcept;"},
{"lineNum":"   59","line":"bool operator> (thread::id x, thread::id y) noexcept;"},
{"lineNum":"   60","line":"bool operator>=(thread::id x, thread::id y) noexcept;"},
{"lineNum":"   61","line":""},
{"lineNum":"   62","line":"template<class charT, class traits>"},
{"lineNum":"   63","line":"basic_ostream<charT, traits>&"},
{"lineNum":"   64","line":"operator<<(basic_ostream<charT, traits>& out, thread::id id);"},
{"lineNum":"   65","line":""},
{"lineNum":"   66","line":"namespace this_thread"},
{"lineNum":"   67","line":"{"},
{"lineNum":"   68","line":""},
{"lineNum":"   69","line":"thread::id get_id() noexcept;"},
{"lineNum":"   70","line":""},
{"lineNum":"   71","line":"void yield() noexcept;"},
{"lineNum":"   72","line":""},
{"lineNum":"   73","line":"template <class Clock, class Duration>"},
{"lineNum":"   74","line":"void sleep_until(const chrono::time_point<Clock, Duration>& abs_time);"},
{"lineNum":"   75","line":""},
{"lineNum":"   76","line":"template <class Rep, class Period>"},
{"lineNum":"   77","line":"void sleep_for(const chrono::duration<Rep, Period>& rel_time);"},
{"lineNum":"   78","line":""},
{"lineNum":"   79","line":"}  // this_thread"},
{"lineNum":"   80","line":""},
{"lineNum":"   81","line":"}  // std"},
{"lineNum":"   82","line":""},
{"lineNum":"   83","line":"*/"},
{"lineNum":"   84","line":""},
{"lineNum":"   85","line":"#include <__config>"},
{"lineNum":"   86","line":"#include <__debug>"},
{"lineNum":"   87","line":"#include <__functional_base>"},
{"lineNum":"   88","line":"#include <__mutex_base>"},
{"lineNum":"   89","line":"#include <__threading_support>"},
{"lineNum":"   90","line":"#include <__utility/__decay_copy.h>"},
{"lineNum":"   91","line":"#include <__utility/forward.h>"},
{"lineNum":"   92","line":"#include <chrono>"},
{"lineNum":"   93","line":"#include <cstddef>"},
{"lineNum":"   94","line":"#include <functional>"},
{"lineNum":"   95","line":"#include <iosfwd>"},
{"lineNum":"   96","line":"#include <memory>"},
{"lineNum":"   97","line":"#include <system_error>"},
{"lineNum":"   98","line":"#include <tuple>"},
{"lineNum":"   99","line":"#include <type_traits>"},
{"lineNum":"  100","line":""},
{"lineNum":"  101","line":"#if !defined(_LIBCPP_HAS_NO_PRAGMA_SYSTEM_HEADER)"},
{"lineNum":"  102","line":"#pragma GCC system_header"},
{"lineNum":"  103","line":"#endif"},
{"lineNum":"  104","line":""},
{"lineNum":"  105","line":"_LIBCPP_PUSH_MACROS"},
{"lineNum":"  106","line":"#include <__undef_macros>"},
{"lineNum":"  107","line":""},
{"lineNum":"  108","line":"#ifdef _LIBCPP_HAS_NO_THREADS"},
{"lineNum":"  109","line":"#error <thread> is not supported on this single threaded system"},
{"lineNum":"  110","line":"#else // !_LIBCPP_HAS_NO_THREADS"},
{"lineNum":"  111","line":""},
{"lineNum":"  112","line":"_LIBCPP_BEGIN_NAMESPACE_STD"},
{"lineNum":"  113","line":""},
{"lineNum":"  114","line":"template <class _Tp> class __thread_specific_ptr;"},
{"lineNum":"  115","line":"class _LIBCPP_TYPE_VIS __thread_struct;"},
{"lineNum":"  116","line":"class _LIBCPP_HIDDEN __thread_struct_imp;"},
{"lineNum":"  117","line":"class __assoc_sub_state;"},
{"lineNum":"  118","line":""},
{"lineNum":"  119","line":"_LIBCPP_FUNC_VIS __thread_specific_ptr<__thread_struct>& __thread_local_data();"},
{"lineNum":"  120","line":""},
{"lineNum":"  121","line":"class _LIBCPP_TYPE_VIS __thread_struct"},
{"lineNum":"  122","line":"{"},
{"lineNum":"  123","line":"    __thread_struct_imp* __p_;"},
{"lineNum":"  124","line":""},
{"lineNum":"  125","line":"    __thread_struct(const __thread_struct&);"},
{"lineNum":"  126","line":"    __thread_struct& operator=(const __thread_struct&);"},
{"lineNum":"  127","line":"public:"},
{"lineNum":"  128","line":"    __thread_struct();"},
{"lineNum":"  129","line":"    ~__thread_struct();"},
{"lineNum":"  130","line":""},
{"lineNum":"  131","line":"    void notify_all_at_thread_exit(condition_variable*, mutex*);"},
{"lineNum":"  132","line":"    void __make_ready_at_thread_exit(__assoc_sub_state*);"},
{"lineNum":"  133","line":"};"},
{"lineNum":"  134","line":""},
{"lineNum":"  135","line":"template <class _Tp>"},
{"lineNum":"  136","line":"class __thread_specific_ptr"},
{"lineNum":"  137","line":"{"},
{"lineNum":"  138","line":"    __libcpp_tls_key __key_;"},
{"lineNum":"  139","line":""},
{"lineNum":"  140","line":"     // Only __thread_local_data() may construct a __thread_specific_ptr"},
{"lineNum":"  141","line":"     // and only with _Tp == __thread_struct."},
{"lineNum":"  142","line":"    static_assert((is_same<_Tp, __thread_struct>::value), \"\");"},
{"lineNum":"  143","line":"    __thread_specific_ptr();"},
{"lineNum":"  144","line":"    friend _LIBCPP_FUNC_VIS __thread_specific_ptr<__thread_struct>& __thread_local_data();"},
{"lineNum":"  145","line":""},
{"lineNum":"  146","line":"    __thread_specific_ptr(const __thread_specific_ptr&);"},
{"lineNum":"  147","line":"    __thread_specific_ptr& operator=(const __thread_specific_ptr&);"},
{"lineNum":"  148","line":""},
{"lineNum":"  149","line":"    _LIBCPP_HIDDEN static void _LIBCPP_TLS_DESTRUCTOR_CC __at_thread_exit(void*);"},
{"lineNum":"  150","line":""},
{"lineNum":"  151","line":"public:"},
{"lineNum":"  152","line":"    typedef _Tp* pointer;"},
{"lineNum":"  153","line":""},
{"lineNum":"  154","line":"    ~__thread_specific_ptr();"},
{"lineNum":"  155","line":""},
{"lineNum":"  156","line":"    _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  157","line":"    pointer get() const {return static_cast<_Tp*>(__libcpp_tls_get(__key_));}"},
{"lineNum":"  158","line":"    _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  159","line":"    pointer operator*() const {return *get();}"},
{"lineNum":"  160","line":"    _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  161","line":"    pointer operator->() const {return get();}"},
{"lineNum":"  162","line":"    void set_pointer(pointer __p);"},
{"lineNum":"  163","line":"};"},
{"lineNum":"  164","line":""},
{"lineNum":"  165","line":"template <class _Tp>"},
{"lineNum":"  166","line":"void _LIBCPP_TLS_DESTRUCTOR_CC"},
{"lineNum":"  167","line":"__thread_specific_ptr<_Tp>::__at_thread_exit(void* __p)"},
{"lineNum":"  168","line":"{"},
{"lineNum":"  169","line":"    delete static_cast<pointer>(__p);"},
{"lineNum":"  170","line":"}"},
{"lineNum":"  171","line":""},
{"lineNum":"  172","line":"template <class _Tp>"},
{"lineNum":"  173","line":"__thread_specific_ptr<_Tp>::__thread_specific_ptr()"},
{"lineNum":"  174","line":"{"},
{"lineNum":"  175","line":"  int __ec ="},
{"lineNum":"  176","line":"      __libcpp_tls_create(&__key_, &__thread_specific_ptr::__at_thread_exit);"},
{"lineNum":"  177","line":"  if (__ec)"},
{"lineNum":"  178","line":"    __throw_system_error(__ec, \"__thread_specific_ptr construction failed\");"},
{"lineNum":"  179","line":"}"},
{"lineNum":"  180","line":""},
{"lineNum":"  181","line":"template <class _Tp>"},
{"lineNum":"  182","line":"__thread_specific_ptr<_Tp>::~__thread_specific_ptr()"},
{"lineNum":"  183","line":"{"},
{"lineNum":"  184","line":"    // __thread_specific_ptr is only created with a static storage duration"},
{"lineNum":"  185","line":"    // so this destructor is only invoked during program termination. Invoking"},
{"lineNum":"  186","line":"    // pthread_key_delete(__key_) may prevent other threads from deleting their"},
{"lineNum":"  187","line":"    // thread local data. For this reason we leak the key."},
{"lineNum":"  188","line":"}"},
{"lineNum":"  189","line":""},
{"lineNum":"  190","line":"template <class _Tp>"},
{"lineNum":"  191","line":"void"},
{"lineNum":"  192","line":"__thread_specific_ptr<_Tp>::set_pointer(pointer __p)"},
{"lineNum":"  193","line":"{"},
{"lineNum":"  194","line":"    _LIBCPP_ASSERT(get() == nullptr,"},
{"lineNum":"  195","line":"                   \"Attempting to overwrite thread local data\");"},
{"lineNum":"  196","line":"    __libcpp_tls_set(__key_, __p);"},
{"lineNum":"  197","line":"}"},
{"lineNum":"  198","line":""},
{"lineNum":"  199","line":"template<>"},
{"lineNum":"  200","line":"struct _LIBCPP_TEMPLATE_VIS hash<__thread_id>"},
{"lineNum":"  201","line":"    : public unary_function<__thread_id, size_t>"},
{"lineNum":"  202","line":"{"},
{"lineNum":"  203","line":"    _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  204","line":"    size_t operator()(__thread_id __v) const _NOEXCEPT"},
{"lineNum":"  205","line":"    {"},
{"lineNum":"  206","line":"        return hash<__libcpp_thread_id>()(__v.__id_);"},
{"lineNum":"  207","line":"    }"},
{"lineNum":"  208","line":"};"},
{"lineNum":"  209","line":""},
{"lineNum":"  210","line":"template<class _CharT, class _Traits>"},
{"lineNum":"  211","line":"_LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  212","line":"basic_ostream<_CharT, _Traits>&"},
{"lineNum":"  213","line":"operator<<(basic_ostream<_CharT, _Traits>& __os, __thread_id __id)"},
{"lineNum":"  214","line":"{return __os << __id.__id_;}"},
{"lineNum":"  215","line":""},
{"lineNum":"  216","line":"class _LIBCPP_TYPE_VIS thread"},
{"lineNum":"  217","line":"{"},
{"lineNum":"  218","line":"    __libcpp_thread_t __t_;"},
{"lineNum":"  219","line":""},
{"lineNum":"  220","line":"    thread(const thread&);"},
{"lineNum":"  221","line":"    thread& operator=(const thread&);"},
{"lineNum":"  222","line":"public:"},
{"lineNum":"  223","line":"    typedef __thread_id id;"},
{"lineNum":"  224","line":"    typedef __libcpp_thread_t native_handle_type;"},
{"lineNum":"  225","line":""},
{"lineNum":"  226","line":"    _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  227","line":"    thread() _NOEXCEPT : __t_(_LIBCPP_NULL_THREAD) {}"},
{"lineNum":"  228","line":"#ifndef _LIBCPP_CXX03_LANG"},
{"lineNum":"  229","line":"    template <class _Fp, class ..._Args,"},
{"lineNum":"  230","line":"              class = typename enable_if"},
{"lineNum":"  231","line":"              <"},
{"lineNum":"  232","line":"                   !is_same<typename __uncvref<_Fp>::type, thread>::value"},
{"lineNum":"  233","line":"              >::type"},
{"lineNum":"  234","line":"             >"},
{"lineNum":"  235","line":"        _LIBCPP_METHOD_TEMPLATE_IMPLICIT_INSTANTIATION_VIS"},
{"lineNum":"  236","line":"        explicit thread(_Fp&& __f, _Args&&... __args);"},
{"lineNum":"  237","line":"#else  // _LIBCPP_CXX03_LANG"},
{"lineNum":"  238","line":"    template <class _Fp>"},
{"lineNum":"  239","line":"    _LIBCPP_METHOD_TEMPLATE_IMPLICIT_INSTANTIATION_VIS"},
{"lineNum":"  240","line":"    explicit thread(_Fp __f);"},
{"lineNum":"  241","line":"#endif"},
{"lineNum":"  242","line":"    ~thread();"},
{"lineNum":"  243","line":""},
{"lineNum":"  244","line":"    _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  245","line":"    thread(thread&& __t) _NOEXCEPT : __t_(__t.__t_) {"},
{"lineNum":"  246","line":"        __t.__t_ = _LIBCPP_NULL_THREAD;"},
{"lineNum":"  247","line":"    }"},
{"lineNum":"  248","line":""},
{"lineNum":"  249","line":"    _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  250","line":"    thread& operator=(thread&& __t) _NOEXCEPT {"},
{"lineNum":"  251","line":"        if (!__libcpp_thread_isnull(&__t_))"},
{"lineNum":"  252","line":"            terminate();"},
{"lineNum":"  253","line":"        __t_ = __t.__t_;"},
{"lineNum":"  254","line":"        __t.__t_ = _LIBCPP_NULL_THREAD;"},
{"lineNum":"  255","line":"        return *this;"},
{"lineNum":"  256","line":"    }"},
{"lineNum":"  257","line":""},
{"lineNum":"  258","line":"    _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  259","line":"    void swap(thread& __t) _NOEXCEPT {_VSTD::swap(__t_, __t.__t_);}"},
{"lineNum":"  260","line":""},
{"lineNum":"  261","line":"    _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  262","line":"    bool joinable() const _NOEXCEPT {return !__libcpp_thread_isnull(&__t_);}"},
{"lineNum":"  263","line":"    void join();"},
{"lineNum":"  264","line":"    void detach();"},
{"lineNum":"  265","line":"    _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  266","line":"    id get_id() const _NOEXCEPT {return __libcpp_thread_get_id(&__t_);}"},
{"lineNum":"  267","line":"    _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  268","line":"    native_handle_type native_handle() _NOEXCEPT {return __t_;}"},
{"lineNum":"  269","line":""},
{"lineNum":"  270","line":"    static unsigned hardware_concurrency() _NOEXCEPT;"},
{"lineNum":"  271","line":"};"},
{"lineNum":"  272","line":""},
{"lineNum":"  273","line":"#ifndef _LIBCPP_CXX03_LANG"},
{"lineNum":"  274","line":""},
{"lineNum":"  275","line":"template <class _TSp, class _Fp, class ..._Args, size_t ..._Indices>"},
{"lineNum":"  276","line":"inline _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  277","line":"void"},
{"lineNum":"  278","line":"__thread_execute(tuple<_TSp, _Fp, _Args...>& __t, __tuple_indices<_Indices...>)"},
{"lineNum":"  279","line":"{"},
{"lineNum":"  280","line":"    _VSTD::__invoke(_VSTD::move(_VSTD::get<1>(__t)), _VSTD::move(_VSTD::get<_Indices>(__t))...);"},
{"lineNum":"  281","line":"}"},
{"lineNum":"  282","line":""},
{"lineNum":"  283","line":"template <class _Fp>"},
{"lineNum":"  284","line":"_LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  285","line":"void* __thread_proxy(void* __vp)"},
{"lineNum":"  286","line":"{"},
{"lineNum":"  287","line":"    // _Fp = tuple< unique_ptr<__thread_struct>, Functor, Args...>"},
{"lineNum":"  288","line":"    unique_ptr<_Fp> __p(static_cast<_Fp*>(__vp));"},
{"lineNum":"  289","line":"    __thread_local_data().set_pointer(_VSTD::get<0>(*__p.get()).release());"},
{"lineNum":"  290","line":"    typedef typename __make_tuple_indices<tuple_size<_Fp>::value, 2>::type _Index;"},
{"lineNum":"  291","line":"    _VSTD::__thread_execute(*__p.get(), _Index());"},
{"lineNum":"  292","line":"    return nullptr;"},
{"lineNum":"  293","line":"}"},
{"lineNum":"  294","line":""},
{"lineNum":"  295","line":"template <class _Fp, class ..._Args,"},
{"lineNum":"  296","line":"          class"},
{"lineNum":"  297","line":"         >"},
{"lineNum":"  298","line":"thread::thread(_Fp&& __f, _Args&&... __args)"},
{"lineNum":"  299","line":"{"},
{"lineNum":"  300","line":"    typedef unique_ptr<__thread_struct> _TSPtr;"},
{"lineNum":"  301","line":"    _TSPtr __tsp(new __thread_struct);"},
{"lineNum":"  302","line":"    typedef tuple<_TSPtr, typename decay<_Fp>::type, typename decay<_Args>::type...> _Gp;"},
{"lineNum":"  303","line":"    unique_ptr<_Gp> __p("},
{"lineNum":"  304","line":"            new _Gp(_VSTD::move(__tsp),"},
{"lineNum":"  305","line":"                    _VSTD::__decay_copy(_VSTD::forward<_Fp>(__f)),"},
{"lineNum":"  306","line":"                    _VSTD::__decay_copy(_VSTD::forward<_Args>(__args))...));"},
{"lineNum":"  307","line":"    int __ec = _VSTD::__libcpp_thread_create(&__t_, &__thread_proxy<_Gp>, __p.get());"},
{"lineNum":"  308","line":"    if (__ec == 0)"},
{"lineNum":"  309","line":"        __p.release();"},
{"lineNum":"  310","line":"    else"},
{"lineNum":"  311","line":"        __throw_system_error(__ec, \"thread constructor failed\");"},
{"lineNum":"  312","line":"}"},
{"lineNum":"  313","line":""},
{"lineNum":"  314","line":"#else  // _LIBCPP_CXX03_LANG"},
{"lineNum":"  315","line":""},
{"lineNum":"  316","line":"template <class _Fp>"},
{"lineNum":"  317","line":"struct __thread_invoke_pair {"},
{"lineNum":"  318","line":"    // This type is used to pass memory for thread local storage and a functor"},
{"lineNum":"  319","line":"    // to a newly created thread because std::pair doesn\'t work with"},
{"lineNum":"  320","line":"    // std::unique_ptr in C++03."},
{"lineNum":"  321","line":"    __thread_invoke_pair(_Fp& __f) : __tsp_(new __thread_struct), __fn_(__f) {}"},
{"lineNum":"  322","line":"    unique_ptr<__thread_struct> __tsp_;"},
{"lineNum":"  323","line":"    _Fp __fn_;"},
{"lineNum":"  324","line":"};"},
{"lineNum":"  325","line":""},
{"lineNum":"  326","line":"template <class _Fp>"},
{"lineNum":"  327","line":"void* __thread_proxy_cxx03(void* __vp)"},
{"lineNum":"  328","line":"{"},
{"lineNum":"  329","line":"    unique_ptr<_Fp> __p(static_cast<_Fp*>(__vp));"},
{"lineNum":"  330","line":"    __thread_local_data().set_pointer(__p->__tsp_.release());"},
{"lineNum":"  331","line":"    (__p->__fn_)();"},
{"lineNum":"  332","line":"    return nullptr;"},
{"lineNum":"  333","line":"}"},
{"lineNum":"  334","line":""},
{"lineNum":"  335","line":"template <class _Fp>"},
{"lineNum":"  336","line":"thread::thread(_Fp __f)"},
{"lineNum":"  337","line":"{"},
{"lineNum":"  338","line":""},
{"lineNum":"  339","line":"    typedef __thread_invoke_pair<_Fp> _InvokePair;"},
{"lineNum":"  340","line":"    typedef unique_ptr<_InvokePair> _PairPtr;"},
{"lineNum":"  341","line":"    _PairPtr __pp(new _InvokePair(__f));"},
{"lineNum":"  342","line":"    int __ec = _VSTD::__libcpp_thread_create(&__t_, &__thread_proxy_cxx03<_InvokePair>, __pp.get());"},
{"lineNum":"  343","line":"    if (__ec == 0)"},
{"lineNum":"  344","line":"        __pp.release();"},
{"lineNum":"  345","line":"    else"},
{"lineNum":"  346","line":"        __throw_system_error(__ec, \"thread constructor failed\");"},
{"lineNum":"  347","line":"}"},
{"lineNum":"  348","line":""},
{"lineNum":"  349","line":"#endif // _LIBCPP_CXX03_LANG"},
{"lineNum":"  350","line":""},
{"lineNum":"  351","line":"inline _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  352","line":"void swap(thread& __x, thread& __y) _NOEXCEPT {__x.swap(__y);}"},
{"lineNum":"  353","line":""},
{"lineNum":"  354","line":"namespace this_thread"},
{"lineNum":"  355","line":"{"},
{"lineNum":"  356","line":""},
{"lineNum":"  357","line":"_LIBCPP_FUNC_VIS void sleep_for(const chrono::nanoseconds& __ns);"},
{"lineNum":"  358","line":""},
{"lineNum":"  359","line":"template <class _Rep, class _Period>"},
{"lineNum":"  360","line":"void"},
{"lineNum":"  361","line":"sleep_for(const chrono::duration<_Rep, _Period>& __d)"},
{"lineNum":"  362","line":"{"},
{"lineNum":"  363","line":"    if (__d > chrono::duration<_Rep, _Period>::zero())","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  364","line":"    {"},
{"lineNum":"  365","line":"        // The standard guarantees a 64bit signed integer resolution for nanoseconds,"},
{"lineNum":"  366","line":"        // so use INT64_MAX / 1e9 as cut-off point. Use a constant to avoid <climits>"},
{"lineNum":"  367","line":"        // and issues with long double folding on PowerPC with GCC."},
{"lineNum":"  368","line":"        _LIBCPP_CONSTEXPR chrono::duration<long double> _Max ="},
{"lineNum":"  369","line":"            chrono::duration<long double>(9223372036.0L);"},
{"lineNum":"  370","line":"        chrono::nanoseconds __ns;"},
{"lineNum":"  371","line":"        if (__d < _Max)"},
{"lineNum":"  372","line":"        {"},
{"lineNum":"  373","line":"            __ns = chrono::duration_cast<chrono::nanoseconds>(__d);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  374","line":"            if (__ns < __d)"},
{"lineNum":"  375","line":"                ++__ns;"},
{"lineNum":"  376","line":"        }"},
{"lineNum":"  377","line":"        else"},
{"lineNum":"  378","line":"            __ns = chrono::nanoseconds::max();"},
{"lineNum":"  379","line":"        this_thread::sleep_for(__ns);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  380","line":"    }"},
{"lineNum":"  381","line":"}"},
{"lineNum":"  382","line":""},
{"lineNum":"  383","line":"template <class _Clock, class _Duration>"},
{"lineNum":"  384","line":"void"},
{"lineNum":"  385","line":"sleep_until(const chrono::time_point<_Clock, _Duration>& __t)"},
{"lineNum":"  386","line":"{"},
{"lineNum":"  387","line":"    mutex __mut;"},
{"lineNum":"  388","line":"    condition_variable __cv;"},
{"lineNum":"  389","line":"    unique_lock<mutex> __lk(__mut);"},
{"lineNum":"  390","line":"    while (_Clock::now() < __t)"},
{"lineNum":"  391","line":"        __cv.wait_until(__lk, __t);"},
{"lineNum":"  392","line":"}"},
{"lineNum":"  393","line":""},
{"lineNum":"  394","line":"template <class _Duration>"},
{"lineNum":"  395","line":"inline _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  396","line":"void"},
{"lineNum":"  397","line":"sleep_until(const chrono::time_point<chrono::steady_clock, _Duration>& __t)"},
{"lineNum":"  398","line":"{"},
{"lineNum":"  399","line":"    this_thread::sleep_for(__t - chrono::steady_clock::now());"},
{"lineNum":"  400","line":"}"},
{"lineNum":"  401","line":""},
{"lineNum":"  402","line":"inline _LIBCPP_INLINE_VISIBILITY"},
{"lineNum":"  403","line":"void yield() _NOEXCEPT {__libcpp_thread_yield();}"},
{"lineNum":"  404","line":""},
{"lineNum":"  405","line":"}  // this_thread"},
{"lineNum":"  406","line":""},
{"lineNum":"  407","line":"_LIBCPP_END_NAMESPACE_STD"},
{"lineNum":"  408","line":""},
{"lineNum":"  409","line":"#endif // !_LIBCPP_HAS_NO_THREADS"},
{"lineNum":"  410","line":""},
{"lineNum":"  411","line":"_LIBCPP_POP_MACROS"},
{"lineNum":"  412","line":""},
{"lineNum":"  413","line":"#endif // _LIBCPP_THREAD"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 3, "covered" : 0,};
var merged_data = [];

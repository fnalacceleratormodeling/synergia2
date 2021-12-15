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
{"lineNum":"   44","line":"#ifndef KOKKOS_ATOMIC_VIEW_HPP"},
{"lineNum":"   45","line":"#define KOKKOS_ATOMIC_VIEW_HPP"},
{"lineNum":"   46","line":""},
{"lineNum":"   47","line":"#include <Kokkos_Macros.hpp>"},
{"lineNum":"   48","line":"#include <Kokkos_Atomic.hpp>"},
{"lineNum":"   49","line":""},
{"lineNum":"   50","line":"namespace Kokkos {"},
{"lineNum":"   51","line":"namespace Impl {"},
{"lineNum":"   52","line":""},
{"lineNum":"   53","line":"// The following tag is used to prevent an implicit call of the constructor when"},
{"lineNum":"   54","line":"// trying to assign a literal 0 int ( = 0 );"},
{"lineNum":"   55","line":"struct AtomicViewConstTag {};"},
{"lineNum":"   56","line":""},
{"lineNum":"   57","line":"template <class ViewTraits>"},
{"lineNum":"   58","line":"class AtomicDataElement {"},
{"lineNum":"   59","line":" public:"},
{"lineNum":"   60","line":"  using value_type           = typename ViewTraits::value_type;"},
{"lineNum":"   61","line":"  using const_value_type     = typename ViewTraits::const_value_type;"},
{"lineNum":"   62","line":"  using non_const_value_type = typename ViewTraits::non_const_value_type;"},
{"lineNum":"   63","line":"  volatile value_type* const ptr;"},
{"lineNum":"   64","line":""},
{"lineNum":"   65","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"   66","line":"  AtomicDataElement(value_type* ptr_, AtomicViewConstTag) : ptr(ptr_) {}"},
{"lineNum":"   67","line":""},
{"lineNum":"   68","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"   69","line":"  const_value_type operator=(const_value_type& val) const {"},
{"lineNum":"   70","line":"    *ptr = val;"},
{"lineNum":"   71","line":"    return val;"},
{"lineNum":"   72","line":"  }"},
{"lineNum":"   73","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"   74","line":"  const_value_type operator=(volatile const_value_type& val) const {"},
{"lineNum":"   75","line":"    *ptr = val;"},
{"lineNum":"   76","line":"    return val;"},
{"lineNum":"   77","line":"  }"},
{"lineNum":"   78","line":""},
{"lineNum":"   79","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"   80","line":"  void inc() const { Kokkos::atomic_increment(ptr); }"},
{"lineNum":"   81","line":""},
{"lineNum":"   82","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"   83","line":"  void dec() const { Kokkos::atomic_decrement(ptr); }"},
{"lineNum":"   84","line":""},
{"lineNum":"   85","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"   86","line":"  const_value_type operator++() const {"},
{"lineNum":"   87","line":"    const_value_type tmp ="},
{"lineNum":"   88","line":"        Kokkos::atomic_fetch_add(ptr, non_const_value_type(1));"},
{"lineNum":"   89","line":"    return tmp + 1;"},
{"lineNum":"   90","line":"  }"},
{"lineNum":"   91","line":""},
{"lineNum":"   92","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"   93","line":"  const_value_type operator--() const {"},
{"lineNum":"   94","line":"    const_value_type tmp ="},
{"lineNum":"   95","line":"        Kokkos::atomic_fetch_sub(ptr, non_const_value_type(1));"},
{"lineNum":"   96","line":"    return tmp - 1;"},
{"lineNum":"   97","line":"  }"},
{"lineNum":"   98","line":""},
{"lineNum":"   99","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  100","line":"  const_value_type operator++(int) const {"},
{"lineNum":"  101","line":"    return Kokkos::atomic_fetch_add(ptr, non_const_value_type(1));"},
{"lineNum":"  102","line":"  }"},
{"lineNum":"  103","line":""},
{"lineNum":"  104","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  105","line":"  const_value_type operator--(int) const {"},
{"lineNum":"  106","line":"    return Kokkos::atomic_fetch_sub(ptr, non_const_value_type(1));"},
{"lineNum":"  107","line":"  }"},
{"lineNum":"  108","line":""},
{"lineNum":"  109","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  110","line":"  const_value_type operator+=(const_value_type& val) const {"},
{"lineNum":"  111","line":"    const_value_type tmp = Kokkos::atomic_fetch_add(ptr, val);"},
{"lineNum":"  112","line":"    return tmp + val;"},
{"lineNum":"  113","line":"  }"},
{"lineNum":"  114","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  115","line":"  const_value_type operator+=(volatile const_value_type& val) const {"},
{"lineNum":"  116","line":"    const_value_type tmp = Kokkos::atomic_fetch_add(ptr, val);"},
{"lineNum":"  117","line":"    return tmp + val;"},
{"lineNum":"  118","line":"  }"},
{"lineNum":"  119","line":""},
{"lineNum":"  120","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  121","line":"  const_value_type operator-=(const_value_type& val) const {"},
{"lineNum":"  122","line":"    const_value_type tmp = Kokkos::atomic_fetch_sub(ptr, val);"},
{"lineNum":"  123","line":"    return tmp - val;"},
{"lineNum":"  124","line":"  }"},
{"lineNum":"  125","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  126","line":"  const_value_type operator-=(volatile const_value_type& val) const {"},
{"lineNum":"  127","line":"    const_value_type tmp = Kokkos::atomic_fetch_sub(ptr, val);"},
{"lineNum":"  128","line":"    return tmp - val;"},
{"lineNum":"  129","line":"  }"},
{"lineNum":"  130","line":""},
{"lineNum":"  131","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  132","line":"  const_value_type operator*=(const_value_type& val) const {"},
{"lineNum":"  133","line":"    return Kokkos::atomic_mul_fetch(ptr, val);"},
{"lineNum":"  134","line":"  }"},
{"lineNum":"  135","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  136","line":"  const_value_type operator*=(volatile const_value_type& val) const {"},
{"lineNum":"  137","line":"    return Kokkos::atomic_mul_fetch(ptr, val);"},
{"lineNum":"  138","line":"  }"},
{"lineNum":"  139","line":""},
{"lineNum":"  140","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  141","line":"  const_value_type operator/=(const_value_type& val) const {"},
{"lineNum":"  142","line":"    return Kokkos::atomic_div_fetch(ptr, val);"},
{"lineNum":"  143","line":"  }"},
{"lineNum":"  144","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  145","line":"  const_value_type operator/=(volatile const_value_type& val) const {"},
{"lineNum":"  146","line":"    return Kokkos::atomic_div_fetch(ptr, val);"},
{"lineNum":"  147","line":"  }"},
{"lineNum":"  148","line":""},
{"lineNum":"  149","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  150","line":"  const_value_type operator%=(const_value_type& val) const {"},
{"lineNum":"  151","line":"    return Kokkos::atomic_mod_fetch(ptr, val);"},
{"lineNum":"  152","line":"  }"},
{"lineNum":"  153","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  154","line":"  const_value_type operator%=(volatile const_value_type& val) const {"},
{"lineNum":"  155","line":"    return Kokkos::atomic_mod_fetch(ptr, val);"},
{"lineNum":"  156","line":"  }"},
{"lineNum":"  157","line":""},
{"lineNum":"  158","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  159","line":"  const_value_type operator&=(const_value_type& val) const {"},
{"lineNum":"  160","line":"    return Kokkos::atomic_and_fetch(ptr, val);"},
{"lineNum":"  161","line":"  }"},
{"lineNum":"  162","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  163","line":"  const_value_type operator&=(volatile const_value_type& val) const {"},
{"lineNum":"  164","line":"    return Kokkos::atomic_and_fetch(ptr, val);"},
{"lineNum":"  165","line":"  }"},
{"lineNum":"  166","line":""},
{"lineNum":"  167","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  168","line":"  const_value_type operator^=(const_value_type& val) const {"},
{"lineNum":"  169","line":"    return Kokkos::atomic_xor_fetch(ptr, val);"},
{"lineNum":"  170","line":"  }"},
{"lineNum":"  171","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  172","line":"  const_value_type operator^=(volatile const_value_type& val) const {"},
{"lineNum":"  173","line":"    return Kokkos::atomic_xor_fetch(ptr, val);"},
{"lineNum":"  174","line":"  }"},
{"lineNum":"  175","line":""},
{"lineNum":"  176","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  177","line":"  const_value_type operator|=(const_value_type& val) const {"},
{"lineNum":"  178","line":"    return Kokkos::atomic_or_fetch(ptr, val);"},
{"lineNum":"  179","line":"  }"},
{"lineNum":"  180","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  181","line":"  const_value_type operator|=(volatile const_value_type& val) const {"},
{"lineNum":"  182","line":"    return Kokkos::atomic_or_fetch(ptr, val);"},
{"lineNum":"  183","line":"  }"},
{"lineNum":"  184","line":""},
{"lineNum":"  185","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  186","line":"  const_value_type operator<<=(const_value_type& val) const {"},
{"lineNum":"  187","line":"    return Kokkos::atomic_lshift_fetch(ptr, val);"},
{"lineNum":"  188","line":"  }"},
{"lineNum":"  189","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  190","line":"  const_value_type operator<<=(volatile const_value_type& val) const {"},
{"lineNum":"  191","line":"    return Kokkos::atomic_lshift_fetch(ptr, val);"},
{"lineNum":"  192","line":"  }"},
{"lineNum":"  193","line":""},
{"lineNum":"  194","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  195","line":"  const_value_type operator>>=(const_value_type& val) const {"},
{"lineNum":"  196","line":"    return Kokkos::atomic_rshift_fetch(ptr, val);"},
{"lineNum":"  197","line":"  }"},
{"lineNum":"  198","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  199","line":"  const_value_type operator>>=(volatile const_value_type& val) const {"},
{"lineNum":"  200","line":"    return Kokkos::atomic_rshift_fetch(ptr, val);"},
{"lineNum":"  201","line":"  }"},
{"lineNum":"  202","line":""},
{"lineNum":"  203","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  204","line":"  const_value_type operator+(const_value_type& val) const { return *ptr + val; }"},
{"lineNum":"  205","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  206","line":"  const_value_type operator+(volatile const_value_type& val) const {"},
{"lineNum":"  207","line":"    return *ptr + val;"},
{"lineNum":"  208","line":"  }"},
{"lineNum":"  209","line":""},
{"lineNum":"  210","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  211","line":"  const_value_type operator-(const_value_type& val) const { return *ptr - val; }"},
{"lineNum":"  212","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  213","line":"  const_value_type operator-(volatile const_value_type& val) const {"},
{"lineNum":"  214","line":"    return *ptr - val;"},
{"lineNum":"  215","line":"  }"},
{"lineNum":"  216","line":""},
{"lineNum":"  217","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  218","line":"  const_value_type operator*(const_value_type& val) const { return *ptr * val; }"},
{"lineNum":"  219","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  220","line":"  const_value_type operator*(volatile const_value_type& val) const {"},
{"lineNum":"  221","line":"    return *ptr * val;"},
{"lineNum":"  222","line":"  }"},
{"lineNum":"  223","line":""},
{"lineNum":"  224","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  225","line":"  const_value_type operator/(const_value_type& val) const { return *ptr / val; }"},
{"lineNum":"  226","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  227","line":"  const_value_type operator/(volatile const_value_type& val) const {"},
{"lineNum":"  228","line":"    return *ptr / val;"},
{"lineNum":"  229","line":"  }"},
{"lineNum":"  230","line":""},
{"lineNum":"  231","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  232","line":"  const_value_type operator%(const_value_type& val) const { return *ptr ^ val; }"},
{"lineNum":"  233","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  234","line":"  const_value_type operator%(volatile const_value_type& val) const {"},
{"lineNum":"  235","line":"    return *ptr ^ val;"},
{"lineNum":"  236","line":"  }"},
{"lineNum":"  237","line":""},
{"lineNum":"  238","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  239","line":"  const_value_type operator!() const { return !*ptr; }"},
{"lineNum":"  240","line":""},
{"lineNum":"  241","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  242","line":"  const_value_type operator&&(const_value_type& val) const {"},
{"lineNum":"  243","line":"    return *ptr && val;"},
{"lineNum":"  244","line":"  }"},
{"lineNum":"  245","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  246","line":"  const_value_type operator&&(volatile const_value_type& val) const {"},
{"lineNum":"  247","line":"    return *ptr && val;"},
{"lineNum":"  248","line":"  }"},
{"lineNum":"  249","line":""},
{"lineNum":"  250","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  251","line":"  const_value_type operator||(const_value_type& val) const {"},
{"lineNum":"  252","line":"    return *ptr | val;"},
{"lineNum":"  253","line":"  }"},
{"lineNum":"  254","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  255","line":"  const_value_type operator||(volatile const_value_type& val) const {"},
{"lineNum":"  256","line":"    return *ptr | val;"},
{"lineNum":"  257","line":"  }"},
{"lineNum":"  258","line":""},
{"lineNum":"  259","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  260","line":"  const_value_type operator&(const_value_type& val) const { return *ptr & val; }"},
{"lineNum":"  261","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  262","line":"  const_value_type operator&(volatile const_value_type& val) const {"},
{"lineNum":"  263","line":"    return *ptr & val;"},
{"lineNum":"  264","line":"  }"},
{"lineNum":"  265","line":""},
{"lineNum":"  266","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  267","line":"  const_value_type operator|(const_value_type& val) const { return *ptr | val; }"},
{"lineNum":"  268","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  269","line":"  const_value_type operator|(volatile const_value_type& val) const {"},
{"lineNum":"  270","line":"    return *ptr | val;"},
{"lineNum":"  271","line":"  }"},
{"lineNum":"  272","line":""},
{"lineNum":"  273","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  274","line":"  const_value_type operator^(const_value_type& val) const { return *ptr ^ val; }"},
{"lineNum":"  275","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  276","line":"  const_value_type operator^(volatile const_value_type& val) const {"},
{"lineNum":"  277","line":"    return *ptr ^ val;"},
{"lineNum":"  278","line":"  }"},
{"lineNum":"  279","line":""},
{"lineNum":"  280","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  281","line":"  const_value_type operator~() const { return ~*ptr; }"},
{"lineNum":"  282","line":""},
{"lineNum":"  283","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  284","line":"  const_value_type operator<<(const unsigned int& val) const {"},
{"lineNum":"  285","line":"    return *ptr << val;"},
{"lineNum":"  286","line":"  }"},
{"lineNum":"  287","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  288","line":"  const_value_type operator<<(volatile const unsigned int& val) const {"},
{"lineNum":"  289","line":"    return *ptr << val;"},
{"lineNum":"  290","line":"  }"},
{"lineNum":"  291","line":""},
{"lineNum":"  292","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  293","line":"  const_value_type operator>>(const unsigned int& val) const {"},
{"lineNum":"  294","line":"    return *ptr >> val;"},
{"lineNum":"  295","line":"  }"},
{"lineNum":"  296","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  297","line":"  const_value_type operator>>(volatile const unsigned int& val) const {"},
{"lineNum":"  298","line":"    return *ptr >> val;"},
{"lineNum":"  299","line":"  }"},
{"lineNum":"  300","line":""},
{"lineNum":"  301","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  302","line":"  bool operator==(const AtomicDataElement& val) const { return *ptr == val; }"},
{"lineNum":"  303","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  304","line":"  bool operator==(volatile const AtomicDataElement& val) const {"},
{"lineNum":"  305","line":"    return *ptr == val;"},
{"lineNum":"  306","line":"  }"},
{"lineNum":"  307","line":""},
{"lineNum":"  308","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  309","line":"  bool operator!=(const AtomicDataElement& val) const { return *ptr != val; }"},
{"lineNum":"  310","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  311","line":"  bool operator!=(volatile const AtomicDataElement& val) const {"},
{"lineNum":"  312","line":"    return *ptr != val;"},
{"lineNum":"  313","line":"  }"},
{"lineNum":"  314","line":""},
{"lineNum":"  315","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  316","line":"  bool operator>=(const_value_type& val) const { return *ptr >= val; }"},
{"lineNum":"  317","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  318","line":"  bool operator>=(volatile const_value_type& val) const { return *ptr >= val; }"},
{"lineNum":"  319","line":""},
{"lineNum":"  320","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  321","line":"  bool operator<=(const_value_type& val) const { return *ptr <= val; }"},
{"lineNum":"  322","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  323","line":"  bool operator<=(volatile const_value_type& val) const { return *ptr <= val; }"},
{"lineNum":"  324","line":""},
{"lineNum":"  325","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  326","line":"  bool operator<(const_value_type& val) const { return *ptr < val; }"},
{"lineNum":"  327","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  328","line":"  bool operator<(volatile const_value_type& val) const { return *ptr < val; }"},
{"lineNum":"  329","line":""},
{"lineNum":"  330","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  331","line":"  bool operator>(const_value_type& val) const { return *ptr > val; }"},
{"lineNum":"  332","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  333","line":"  bool operator>(volatile const_value_type& val) const { return *ptr > val; }"},
{"lineNum":"  334","line":""},
{"lineNum":"  335","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  336","line":"  operator const_value_type() const {"},
{"lineNum":"  337","line":"    // return Kokkos::atomic_load(ptr);"},
{"lineNum":"  338","line":"    return *ptr;"},
{"lineNum":"  339","line":"  }"},
{"lineNum":"  340","line":""},
{"lineNum":"  341","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  342","line":"  operator volatile non_const_value_type() volatile const {"},
{"lineNum":"  343","line":"    // return Kokkos::atomic_load(ptr);"},
{"lineNum":"  344","line":"    return *ptr;"},
{"lineNum":"  345","line":"  }"},
{"lineNum":"  346","line":"};"},
{"lineNum":"  347","line":""},
{"lineNum":"  348","line":"template <class ViewTraits>"},
{"lineNum":"  349","line":"class AtomicViewDataHandle {"},
{"lineNum":"  350","line":" public:"},
{"lineNum":"  351","line":"  typename ViewTraits::value_type* ptr;"},
{"lineNum":"  352","line":""},
{"lineNum":"  353","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  354","line":"  AtomicViewDataHandle() : ptr(nullptr) {}"},
{"lineNum":"  355","line":""},
{"lineNum":"  356","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  357","line":"  AtomicViewDataHandle(typename ViewTraits::value_type* ptr_) : ptr(ptr_) {}"},
{"lineNum":"  358","line":""},
{"lineNum":"  359","line":"  template <class iType>"},
{"lineNum":"  360","line":"  KOKKOS_INLINE_FUNCTION AtomicDataElement<ViewTraits> operator[]("},
{"lineNum":"  361","line":"      const iType& i) const {"},
{"lineNum":"  362","line":"    return AtomicDataElement<ViewTraits>(ptr + i, AtomicViewConstTag());","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  363","line":"  }"},
{"lineNum":"  364","line":""},
{"lineNum":"  365","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  366","line":"  operator typename ViewTraits::value_type*() const { return ptr; }"},
{"lineNum":"  367","line":"};"},
{"lineNum":"  368","line":""},
{"lineNum":"  369","line":"template <unsigned Size>"},
{"lineNum":"  370","line":"struct Kokkos_Atomic_is_only_allowed_with_32bit_and_64bit_scalars;"},
{"lineNum":"  371","line":""},
{"lineNum":"  372","line":"template <>"},
{"lineNum":"  373","line":"struct Kokkos_Atomic_is_only_allowed_with_32bit_and_64bit_scalars<4> {"},
{"lineNum":"  374","line":"  using type = int;"},
{"lineNum":"  375","line":"};"},
{"lineNum":"  376","line":""},
{"lineNum":"  377","line":"template <>"},
{"lineNum":"  378","line":"struct Kokkos_Atomic_is_only_allowed_with_32bit_and_64bit_scalars<8> {"},
{"lineNum":"  379","line":"  using type = int64_t;"},
{"lineNum":"  380","line":"};"},
{"lineNum":"  381","line":""},
{"lineNum":"  382","line":"}  // namespace Impl"},
{"lineNum":"  383","line":"}  // namespace Kokkos"},
{"lineNum":"  384","line":""},
{"lineNum":"  385","line":"#endif"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 1, "covered" : 0,};
var merged_data = [];

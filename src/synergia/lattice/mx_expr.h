#ifndef MX_EXPR_H
#define MX_EXPR_H

#include <string>
#include <utility>
#include <vector>
#include <boost/variant.hpp>

namespace synergia
{
  typedef double(*nfunc_t)();
  typedef double(*ufunc_t)(double);
  typedef double(*bfunc_t)(double, double);

  typedef std::pair<std::string, std::string> string_pair_t;

  struct nop_t;
  struct uop_t;
  struct bop_t;

  class MadX;
  class mx_calculator;

  // whether the expr contains a ref string
  class mx_expr_ref_checker;

  typedef boost::variant< double
                        , std::string
                        , string_pair_t
                        , boost::recursive_wrapper<nop_t>
                        , boost::recursive_wrapper<uop_t>
                        , boost::recursive_wrapper<bop_t>
                        > mx_expr;

  typedef std::vector<mx_expr> mx_exprs;

  // evaluate the expression
  double mx_eval(mx_expr const & expr);
  double mx_eval(mx_expr const & expr, MadX const & mx);

  // whether can be evaluated to a number
  bool mx_expr_is_number(mx_expr const& expr);

  // retrieve the ref as a string
  std::string mx_expr_refstr(mx_expr const & expr);
}

struct synergia::nop_t
{
  nop_t(mx_expr const & e)
    : expr(e) { }

  mx_expr expr;
};

struct synergia::uop_t
{
  uop_t(ufunc_t f, mx_expr const & param)
    : param(param), func(f) { }

  mx_expr param;
  ufunc_t func;
};

struct synergia::bop_t
{
  bop_t(bfunc_t f, mx_expr const & lhs, mx_expr const & rhs)
    : func(f), lhs(lhs), rhs(rhs) { }

  bfunc_t func;
  mx_expr lhs;
  mx_expr rhs;
};

class synergia::mx_calculator
  : public boost::static_visitor<double>
{
public:
  mx_calculator() 
    : mx(NULL), def(nan) { }

  explicit 
  mx_calculator(double def)
    : mx(NULL), def(def) { }

  explicit 
  mx_calculator(MadX const & mx) 
    : mx(&mx),  def(nan) { }

  mx_calculator(MadX const & mx, double def) 
    : mx(&mx),  def(def) { }

  double operator()(double val) const;
  double operator()(std::string const & ref) const;
  double operator()(string_pair_t const & ref) const;
  double operator()(nop_t const & n) const;
  double operator()(uop_t const & u) const;
  double operator()(bop_t const & b) const;

public:
  static double nan;

private:
  MadX const * mx;
  double def;
};

class synergia::mx_expr_ref_checker
  : public boost::static_visitor<bool>
{
public:
  bool operator()(double val) const;
  bool operator()(std::string const & ref) const;
  bool operator()(string_pair_t const & ref) const;
  bool operator()(nop_t const & n) const;
  bool operator()(uop_t const & u) const;
  bool operator()(bop_t const & b) const;
};

#endif

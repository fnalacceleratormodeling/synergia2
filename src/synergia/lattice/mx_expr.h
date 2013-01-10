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

  typedef boost::variant< double
                        , std::string
                        , string_pair_t
                        , boost::recursive_wrapper<nop_t>
                        , boost::recursive_wrapper<uop_t>
                        , boost::recursive_wrapper<bop_t>
                        > mx_expr;

  typedef std::vector<mx_expr> mx_exprs;
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
  mx_calculator() : mx(NULL) { }
  mx_calculator(MadX const & mx) : mx(&mx) { }

  double operator()(double val) const;
  double operator()(std::string const & ref) const;
  double operator()(string_pair_t const & ref) const;
  double operator()(nop_t const & n) const;
  double operator()(uop_t const & u) const;
  double operator()(bop_t const & b) const;

private:
  MadX const * mx;
};

#endif

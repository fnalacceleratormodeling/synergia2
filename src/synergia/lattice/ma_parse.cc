
#include "ma_parse.h"

#include <cmath>
#include <limits>
#include <iterator>
#include <fstream>
#include <sstream>

#include <boost/any.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>


namespace ascii = ::boost::spirit::ascii;
namespace phx   = ::boost::phoenix;
namespace qi    = ::boost::spirit::qi;
namespace ql    = ::boost::spirit::qi::labels;

using ascii::char_;
using ascii::no_case;
using ascii::digit;
using ascii::graph;
using ascii::space;

//using phx::ref;
//using phx::bind;

using qi::_1;
using qi::_2;
using qi::_3;
using qi::_val;
using qi::eol;
using qi::lexeme;
using qi::lit;
using qi::bool_;
using qi::int_;
using qi::double_;
using qi::no_skip;
using qi::raw;
using qi::skip;
using qi::locals;

using boost::spirit::qi::real_parser;
using boost::spirit::qi::real_policies;

//using namespace qi::labels;

using namespace std;


// ------------------------------------------------------------------

namespace 
{ 

  struct lazy_pow_
  {
    template <typename X, typename Y>
    struct result { typedef X type; };

    template <typename X, typename Y>
    X operator()(X x, Y y) const
    {
      return std::pow(x, y);
    }
  };

  struct lazy_ufunc_
  {
    template <typename F, typename A1>
    struct result { typedef A1 type; };

    template <typename F, typename A1>
    A1 operator()(F f, A1 a1) const
    {
      return f(a1);
    }
  };

  struct lazy_bfunc_
  {
    template <typename F, typename A1, typename A2>
    struct result { typedef A1 type; };

    template <typename F, typename A1, typename A2>
    A1 operator()(F f, A1 a1, A2 a2) const
    {
      return f(a1, a2);
    }
  };

} 


namespace synergia
{

  template <typename FPT, typename Iterator>
    struct grammar;

  template <typename FPT, typename Iterator>
    bool parse( Iterator &iter,
                Iterator end,
                const grammar<FPT,Iterator> &g,
                FPT &result);

  template <typename Iterator>
    struct mad8_parser;
}


template <typename FPT, typename Iterator>
struct synergia::grammar
  : qi::grammar< Iterator, FPT(), ascii::space_type >
{
  // real policy to allow 'd' and 'D' as exponent
  struct ts_real_policies 
      : boost::spirit::qi::real_policies<FPT>
  {
    static bool 
      parse_exp(Iterator & first, Iterator const &)
    {
      if( string("eEdD").find(*first) != string::npos )
      {
        ++first; return true;
      }

      return false;
    }
  };
        

  // symbol table for constants like "pi"
  struct constant_
    : qi::symbols< typename std::iterator_traits<Iterator>::value_type, 
                   FPT >
  {
    constant_()
    {
      this->add
          ("pi", boost::math::constants::pi<FPT>()  )
      ;
    }
  } constant;

  // symbol table for unary functions like "abs"
  struct ufunc_
    : qi::symbols< typename std::iterator_traits<Iterator>::value_type, 
                   FPT (*)(FPT) >
  {
    ufunc_()
    {
      this->add
          ("abs"   , (FPT (*)(FPT)) std::abs  )
          ("acos"  , (FPT (*)(FPT)) std::acos )
          ("asin"  , (FPT (*)(FPT)) std::asin )
          ("atan"  , (FPT (*)(FPT)) std::atan )
          ("ceil"  , (FPT (*)(FPT)) std::ceil )
          ("cos"   , (FPT (*)(FPT)) std::cos  )
          ("cosh"  , (FPT (*)(FPT)) std::cosh )
          ("exp"   , (FPT (*)(FPT)) std::exp  )
          ("floor" , (FPT (*)(FPT)) std::floor)
          ("log"   , (FPT (*)(FPT)) std::log  )
          ("log10" , (FPT (*)(FPT)) std::log10)
          ("sin"   , (FPT (*)(FPT)) std::sin  )
          ("sinh"  , (FPT (*)(FPT)) std::sinh )
          ("sqrt"  , (FPT (*)(FPT)) std::sqrt )
          ("tan"   , (FPT (*)(FPT)) std::tan  )
          ("tanh"  , (FPT (*)(FPT)) std::tanh )
      ;
    }
  } ufunc;

  // symbol table for binary functions like "pow"
  struct bfunc_
    : qi::symbols< typename std::iterator_traits<Iterator>::value_type, 
                   FPT (*)(FPT, FPT) >
  {
    bfunc_()
    {
      this->add
          ("pow"  , (FPT (*)(FPT, FPT)) std::pow  )
          ("atan2", (FPT (*)(FPT, FPT)) std::atan2)
      ;
    }
  } bfunc;

  void get_ref( FPT & r, string const & key )
  { r = m8.variable_as_number(key); }

  void get_attr_ref( FPT & r, string const & label, string const & key )
  { r = m8.label(label).attribute_as_number(key); }

  qi::rule< Iterator, FPT(), ascii::space_type > expression;
  qi::rule< Iterator, FPT(), ascii::space_type > term;
  qi::rule< Iterator, FPT(), ascii::space_type > factor;
  qi::rule< Iterator, FPT(), ascii::space_type > primary;
  qi::rule< Iterator, string(), ascii::space_type > name;
  qi::real_parser< FPT, ts_real_policies >     real;

  Mad8 & m8;

  grammar(Mad8 & m8) : grammar::base_type(expression), m8(m8)
  {
    boost::phoenix::function<lazy_pow_>   lazy_pow;
    boost::phoenix::function<lazy_ufunc_> lazy_ufunc;
    boost::phoenix::function<lazy_bfunc_> lazy_bfunc;

    expression =
        term                   [_val =  _1]
        >> *(  ('+' >> term    [_val += _1])
            |  ('-' >> term    [_val -= _1])
            )
        ;

    term =
        factor                 [_val =  _1]
        >> *(  ('*' >> factor  [_val *= _1])
            |  ('/' >> factor  [_val /= _1])
            )
        ;

    factor =
        primary                [_val =  _1]
        >> *(  ("^" >> factor  [_val = lazy_pow(_val, _1)])
            )
        ;

    primary =
        real                    [_val =  _1]
        |   '(' >> expression   [_val =  _1] >> ')'
        |   ('-' >> primary     [_val = -_1])
        |   ('+' >> primary     [_val =  _1])
        |   no_case[constant]   [_val =  _1]
        |   (no_case[ufunc] >> '(' >> expression >> ')')
                                [_val = lazy_ufunc(_1, _2)]
        |   (no_case[bfunc] >> '(' >> expression >> ','
                                   >> expression >> ')')
                                [_val = lazy_bfunc(_1, _2, _3)]
        |   (name >> '[' >> name >> ']')
                                [phx::bind(&grammar::get_attr_ref, this, _val, _1, _2)] 
        |   name                [phx::bind(&grammar::get_ref, this, _val, _1)]
        ;

    name = char_("a-zA-Z_") >> *char_(".a-zA-Z_0-9");

  }
};

template <typename FPT, typename Iterator>
bool synergia::parse( Iterator &iter,
                      Iterator end,
                      const grammar<FPT,Iterator> &g,
                      FPT &result)
{
  return boost::spirit::qi::phrase_parse(
              iter, end, g, boost::spirit::ascii::space, result);
}


// semantic actions for mad8 parser
namespace synergia 
{
  namespace 
  {
    void set_cmd_name(Mad8_command & cmd, string const & name)
    { cmd.set_name(name); }

    void ins_cmd_attr_num(Mad8_command & cmd, string const & name, double v)
    { cmd.insert_attribute(name, v); }

    void ins_cmd_attr_str(Mad8_command & cmd, string const & name, string const & v)
    { cmd.insert_attribute(name, v); }

    void ins_cmd_attr_nul(Mad8_command & cmd, string const & name)
    { cmd.insert_attribute(name, string()); }

    void ins_line_op(Mad8_line & line, string const & op)
    { line.insert_operator(op); }

    void ins_line_ele(Mad8_line & line, string const & ele)
    { line.insert_element(ele); }

    void ins_line_line(Mad8_line & line, Mad8_line const & subline)
    { line.insert_subline(subline); }
  }
}

// mad8 parser

template <typename Iterator>
struct synergia::mad8_parser
  : qi::grammar< Iterator, ascii::space_type >
{
  typedef ascii::space_type ws_t;

  // keywords
  struct keywords_
    : boost::spirit::qi::symbols< typename std::iterator_traits<Iterator>::value_type,
                                  string > 
  {
    keywords_() 
    {
      this->add
           ("proton"     , "proton"     )
           ("electron"   , "electron"   )
           ("positron"   , "positron"   )
           ("anti-proton", "anti-proton")
      ;
    }
  } keywords;

  // rules
  qi::rule<Iterator,                 ws_t> statement;
  qi::rule<Iterator,                 ws_t> var_num;
  qi::rule<Iterator,                 ws_t> var_str;
  qi::rule<Iterator,                 ws_t> command;
  qi::rule<Iterator, Mad8_command(), ws_t> cmd;
  qi::rule<Iterator,                 ws_t> line;
  qi::rule<Iterator, Mad8_line(),    ws_t> line_seq;
  qi::rule<Iterator, string(),       ws_t> ops;
  qi::rule<Iterator, string(),       ws_t> identifier;
  qi::rule<Iterator, string(),       ws_t> label;
  qi::rule<Iterator, string(),       ws_t> name;
  qi::rule<Iterator, string(),       ws_t> str;

  grammar<double, Iterator> expr;
  Mad8 & m8;

  void ins_var_num(string const & name, double value)
  { m8.insert_variable(name, value); }

  void ins_var_str(string const & name, string const & value)
  { m8.insert_variable(name, value); }

  void ins_label(string const & name, Mad8_command const & cmd)
  { m8.insert_label(name, cmd); }

  void ins_line(string const & name, Mad8_line const & line)
  { m8.insert_line(name, line); }

  void ins_cmd(Mad8_command const & cmd)
  { m8.insert_command(cmd); }

  // mad8_statement_parser
  mad8_parser(Mad8 & m8) : mad8_parser::base_type( statement ), expr(m8), m8(m8)
  {
    statement = var_num | var_str | line | label | command;

    var_num   = (name >> -lit(':') >> '=' >> expr)
                            [phx::bind(&mad8_parser::ins_var_num, this, _1, _2)]
              ;

    var_str   = (name >> -lit(':') >> '=' >> str)
                            [phx::bind(&mad8_parser::ins_var_str, this, _1, _2)]
              ;

    line      = (name >> ':' >> no_case["line"] >> '=' >> line_seq)
                            [phx::bind(&mad8_parser::ins_line, this, _1, _2)]
              ;

    label     = (name >> ':' >> cmd)
                            [phx::bind(&mad8_parser::ins_label, this, _1, _2)]
              ;

    command   = (cmd)       [phx::bind(&mad8_parser::ins_cmd, this, _1)]
              ;
                      
    cmd       = name        [phx::bind(&set_cmd_name, _val, _1)]
              >> -( ',' >> 
                       ( ( (name >> '=' >> str)  [phx::bind(&ins_cmd_attr_str, _val, _1, _2)]
                         | (name >> '=' >> no_case[keywords])  
                                                 [phx::bind(&ins_cmd_attr_str, _val, _1, _2)]
                         | (name >> '=' >> expr) [phx::bind(&ins_cmd_attr_num, _val, _1, _2)]
                         | (name)                [phx::bind(&ins_cmd_attr_nul, _val, _1)]
                         )
                         % -lit(',') 
                       )
                  )
              ;

    line_seq  = lit('(') >> ( ( -( ops          [phx::bind(&ins_line_op, _val, _1)] )
                                >> ( name       [phx::bind(&ins_line_ele, _val, _1)]
                                     | line_seq [phx::bind(&ins_line_line, _val, _1)]
                                   )
                              ) % ',' ) >> ')'
              ;

    ops       = raw[(int_ >> '*') | '-'];

    identifier = qi::lexeme[char_("a-zA-Z_") >> *char_(".a-zA-Z_0-9") - keywords];
    name = identifier;

    str = qi::lexeme['"' >> +(char_ - '"') >> '"'];  // quoted string
  }

};


bool synergia::parse_mad8_statement( string const & s, Mad8 & m8 )
{
  typedef string::const_iterator iter_t;

  if( s.empty() ) return true;

  mad8_parser<iter_t> parser(m8);

  iter_t       begin = s.begin();
  iter_t const end   = s.end();

  bool b = qi::phrase_parse( begin, end, parser, space )
         && begin == end;

  return b;
}

bool synergia::parse_mad8( string_t const & str, Mad8 & m8 )
{
  stringstream ss(str);

  int nl  = 0;
  int nst = 0;

  bool result;

  string line;
  string statements;

  size_t npos = string_t::npos;

  while( ss.good() )
  {
    getline(ss, line); ++nl;
    //cout << "line " << nl << " ";

    line = line.substr( 0, line.find_first_of('!') );
    //cout << ", '" << line << "'";

    size_t pos = line.find_first_of('&');

    if( pos != npos )
    {
      line = line.substr( 0, pos );
      boost::trim(line);
      statements.append(line);
      //cout << "\n";

      continue;
    }
    else
    {
      boost::trim(line);
      statements.append(line);
      //cout << ", '" << statements << "'";

      vector<string> vs;
      boost::algorithm::split(vs, statements, boost::is_any_of(";"), boost::token_compress_on);

      for(vector<string>::iterator i=vs.begin(); i!=vs.end(); ++i)
      {
        result = synergia::parse_mad8_statement( *i, m8 );
        //cout << ", " << result << "\n";

        if( !result )  { /*cout << *i << "\n";*/ return false; }
      }

      statements.clear();
    }
  }

  return true;
}


bool synergia::parse_mad8_file( string_t const & fname, Mad8 & m8 )
{
  
  ifstream file;
  file.open(fname.c_str());

  if( !file.is_open() )
    throw runtime_error( "Failed to open file " + fname + " for parsing");

  file.seekg(0, std::ios::end);   
  string str;
  str.reserve(file.tellg());
  file.seekg(0, std::ios::beg);

  str.assign((istreambuf_iterator<char>(file)),
              istreambuf_iterator<char>());

  file.close();

  return parse_mad8( str, m8 );
}








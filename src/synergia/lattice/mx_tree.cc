
#include "mx_tree.h"
#include "mx_parse.h"

#include <iostream>
#include <stdexcept>

using namespace synergia;
using namespace std;

using boost::any;
using boost::any_cast;


// helper
namespace synergia
{
  namespace
  {
    // deduct an expression and replace it with a number if all through
    mx_expr simplify( mx_expr const & e, MadX const & mx )
    {
      try
      {
        double r = boost::apply_visitor( mx_calculator(mx), e );
        return mx_expr(r);
      }
      catch(...)
      {
        return e;
      }
    }

    template<typename T>
    void insert_attr(T & t, mx_attr const & attr, MadX const & mx)
    {
      if( attr.type() == MX_ATTR_STRING )
      {
        t.insert_attribute( attr.name(), any_cast<string>(attr.value()) );
      }
      else if( attr.type() == MX_ATTR_PREDEFINED )
      {
        t.insert_attribute( attr.name(), any_cast<mx_keyword>(attr.value()).name );
      }
      else if( attr.type() == MX_ATTR_NUMBER )
      {
        mx_expr e = simplify( any_cast<mx_expr>(attr.value()), mx );
        t.insert_attribute( attr.name(), e );
      }
      else if( attr.type() == MX_ATTR_LAZY_NUMBER )
      {
        t.insert_attribute( attr.name(), any_cast<mx_expr>(attr.value()) );
      }
      else if( attr.type() == MX_ATTR_ARRAY )
      {
        mx_exprs es = any_cast<mx_exprs>(attr.value());
        for( mx_exprs::iterator it = es.begin()
           ; it != es.end(); ++it )
        {
          mx_expr e = simplify( *it, mx );
          *it = e;
        }
        t.insert_attribute( attr.name(), es );
      }
      else if( attr.type() == MX_ATTR_LAZY_ARRAY )
      {
        mx_exprs es = any_cast<mx_exprs>(attr.value());
        t.insert_attribute( attr.name(), es );
      }
    }

  }
}



// attributes
void mx_attr::set_attr(std::string const & name, boost::any const & val)
{
  name_ = name;
  value_ = val;

  if( val.type() == typeid(std::string) )
    type_ = MX_ATTR_STRING;
  else if( val.type() == typeid(mx_expr) )
    type_ = MX_ATTR_NUMBER;
  else if( val.type() == typeid(mx_exprs) )
    type_ = MX_ATTR_ARRAY;
  else if(  val.type() == typeid(mx_keyword)
         && (  any_cast<mx_keyword>(val).tag == MX_KW_PARTICLE
            || any_cast<mx_keyword>(val).tag == MX_KW_MP_TYPE  ) )
    type_ = MX_ATTR_PREDEFINED;
  else
    throw std::runtime_error("Unknown attribute value type for " + name);
}

void mx_attr::set_lazy_attr(std::string const & name, boost::any const & val)
{
  name_ = name;
  value_ = val;

  if( val.type() == typeid(std::string) )
    type_ = MX_ATTR_STRING;      // strings are immediate
  else if( val.type() == typeid(mx_expr) )
    type_ = MX_ATTR_LAZY_NUMBER;
  else if( val.type() == typeid(mx_exprs) )
    type_ = MX_ATTR_LAZY_ARRAY;
  else if(  val.type() == typeid(mx_keyword)
         && (  any_cast<mx_keyword>(val).tag == MX_KW_PARTICLE
            || any_cast<mx_keyword>(val).tag == MX_KW_MP_TYPE  ) )
    type_ = MX_ATTR_PREDEFINED;  // all predefines are immediate
  else
    throw std::runtime_error("Unknown lazy attribute value type for " + name);
}


// command
void mx_command::set_label(string const & label)
{
  label_   = label;
  labeled_ = true;
}

void mx_command::set_keyword(mx_keyword const & keyword)
{
  if( keyword.tag == MX_KW_NONE || keyword.tag == MX_KW_PARTICLE )
    throw std::runtime_error("Invalid keyword type");

  mx_cmd_type t = (keyword.tag == MX_KW_ELEMENT) ? MX_CMD_ELEMENT
                : (keyword.tag == MX_KW_COMMAND) ? MX_CMD_EXECUTABLE : MX_CMD_ELEMENT_REF;

  set_keyword( keyword.name, t );
}

void mx_command::set_keyword(string const & keyword, mx_cmd_type tag)
{
  if( tag!=MX_CMD_ELEMENT && tag!=MX_CMD_EXECUTABLE && tag!=MX_CMD_ELEMENT_REF )
    throw std::runtime_error("Unknown keyword type");

  keyed_   = true;
  keyword_ = keyword;
  type_    = tag;
}

void mx_command::ins_attr(mx_attr const & attr)
{
  attrs_.push_back( attr );
}

bool mx_command::interpret(MadX & mx) const
{
  if( type_ == MX_CMD_VARIABLE )
  {
    mx_attr attr = attrs_[0];
    insert_attr(mx, attr, mx);
  }
  else
  {
    MadX_command_type type =
        (type_ == MX_CMD_ELEMENT) ? ELEMENT
                                  : (type_==MX_CMD_ELEMENT_REF) ? ELEMENT_REF
                                                                : EXECUTABLE ;

    MadX_command cmd;
    cmd.set_name( keyword_, type );
    if ( labeled_ ) {
        cmd.set_label( label_ );
    } else {
        cmd.set_label( "" );
    }
    for( attrs_t::const_iterator it = attrs_.begin()
       ; it != attrs_.end(); ++it )
    {
      mx_attr attr = *it;
      insert_attr(cmd, attr, mx);
    }

    // insert the command to the MadX object
    if( labeled_ ) mx.insert_label(label_, cmd);
    else           mx.insert_command(cmd);
  }

  return true;
}

void mx_command::print() const
{
  if( labeled_ )
    cout << label_ << " : ";

  if( keyed_ )
    cout << keyword_ << ", ";

  for( attrs_t::const_iterator it = attrs_.begin();
       it != attrs_.end(); ++it )
  {
    cout << it->name() << " = " << "xxx, ";
  }

  cout << "\n";
}

// if_block
bool mx_if_block::evaluate_logic(MadX & mx) const
{
  // TODO: parse logic

  return true;
}

bool mx_if_block::interpret_block(MadX & mx) const
{
  return block.interpret(mx);
}

void mx_if_block::print_logic() const
{
  cout << logic_expr;
}

void mx_if_block::print_block() const
{
  block.print();
}

// if-elseif-else
void mx_if::assign_if(string const & logic, mx_tree const & block)
{
  if_ = mx_if_block(logic, block);
}

void mx_if::assign_elseif(string const & logic, mx_tree const & block)
{
  elseif_.push_back(mx_if_block(logic, block));
}

void mx_if::assign_else(mx_tree const & block)
{
  else_ = mx_if_block("1", block);
}

bool mx_if::interpret(MadX & mx) const
{
  if( !if_.valid() )
    throw runtime_error("mx_if::interpret() Invalid if command");

  if( if_.evaluate_logic(mx) )
  {
    return if_.interpret_block(mx);
  }
  else if( !elseif_.empty() )
  {
    for( mx_if_block_v::const_iterator it = elseif_.begin();
         it != elseif_.end();
         ++it )
    {
      if( it->evaluate_logic(mx) )
      {
        return it->interpret_block(mx);
      }
    }
  }
  else if( else_.valid() )
  {
    return else_.interpret_block(mx);
  }

  return true;
}

void mx_if::print() const
{
  if( !if_.valid() )
    throw runtime_error("mx_if::print() Invalid if command");

  cout << "if (";
  if_.print_logic();
  cout << ")\n{\n";
  if_.print_block();
  cout << "}\n";

  for( mx_if_block_v::const_iterator it = elseif_.begin();
       it != elseif_.end(); ++it )
  {
    cout << "elseif (";
    it->print_logic();
    cout << ")\n{\n";
    it->print_block();
    cout << "}\n";
  }

  if( else_.valid() )
  {
    cout << "else\n{\n";
    else_.print_block();
    cout << "}\n";
  }
}

// while
void mx_while::assign(string const & logic, mx_tree const & block)
{
  while_ = mx_if_block(logic, block);
}

bool mx_while::interpret(MadX & mx) const
{
  return true;
}

void mx_while::print() const
{
  if( !while_.valid() )
    throw runtime_error("mx_while::print() Invalid while command");

  cout << "while (";
  while_.print_logic();
  cout << ")\n{\n";
  while_.print_block();
  cout << "}\n";
}

// statement
void mx_statement::assign(mx_command const & st)
{
  value = any(st);
  type = MX_COMMAND;
}

void mx_statement::assign(mx_if const & st)
{
  value = any(st);
  type = MX_IF;
}

void mx_statement::assign(mx_while const & st)
{
  value = any(st);
  type = MX_WHILE;
}

bool mx_statement::interpret(MadX & mx) const
{
  if( type==MX_COMMAND )
    return boost::any_cast<mx_command>(value).interpret(mx);
  else if( type==MX_IF )
    return boost::any_cast<mx_if>(value).interpret(mx);
  else if( type==MX_WHILE )
    return boost::any_cast<mx_while>(value).interpret(mx);
  else
    throw runtime_error("mx_statement::interpret()  Unknown statement type");
}

void mx_statement::print() const
{
  if( type==MX_COMMAND )
    return boost::any_cast<mx_command>(value).print();
  else if( type==MX_IF )
    return boost::any_cast<mx_if>(value).print();
  else if( type==MX_WHILE )
    return boost::any_cast<mx_while>(value).print();
  else
    throw runtime_error("mx_statement::print()  Unknown statement type");
}

// tree
bool mx_tree::interpret(MadX & mx) const
{
  for( mx_statements_t::const_iterator it = statements.begin();
       it != statements.end();
       ++it )
  {
    if( !it->interpret(mx) )
      throw runtime_error("mx_tree::interpret() Failed");
  }
  return true;
}

void mx_tree::print() const
{
  for_each(statements.begin(), statements.end(), detail::print<mx_statement>);
}


















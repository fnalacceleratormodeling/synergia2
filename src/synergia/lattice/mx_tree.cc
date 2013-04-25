
#include "mx_tree.h"
#include "mx_parse.h"

#include <iostream>
#include <stdexcept>
#include <cassert>

#include <boost/algorithm/string.hpp>

#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/four_momentum.h"

using namespace synergia;
using namespace std;

using boost::any;
using boost::any_cast;
using boost::get;


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


void mx_logic::set(mx_expr const & l, logic_op_t o, mx_expr const & r)
{
  lhs = l; rhs = r; op = o; use_preset = false;
}

bool mx_logic::evaluate(MadX const & mx) const
{
  if( use_preset ) return pre;

  assert( op!=NULL );

  double l = boost::apply_visitor( mx_calculator(mx, 0.0), lhs );
  double r = boost::apply_visitor( mx_calculator(mx, 0.0), rhs );

  return op(l, r);
}

// attributes
void mx_attr::set_attr(std::string const & name, boost::any const & val)
{
  name_ = name;
  value_ = val;

  transform(name_.begin(), name_.end(), name_.begin(), ::tolower);

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

  transform(name_.begin(), name_.end(), name_.begin(), ::tolower);

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


void mx_line::interpret(MadX & mx)
{
  MadX_line new_line(mx);
  seq.interpret(mx, new_line, 1);
  mx.insert_line(name, new_line);
}

void mx_line_seq::insert_member(int op, mx_line_member const & member)
{
  members.push_back( make_pair(member, op) );
}

void mx_line_seq::interpret(MadX const & mx, MadX_line & line, int op)
{
  if( op==0 ) return;

  if( op>0 )  // repeat op times
  {
    for(int z=0; z<op; ++z)
    {
      for( mx_line_members::const_iterator it = members.begin()
         ; it!=members.end(); ++it )
      {
        mx_line_member member = it->first;
        int op = it->second;
        member.interpret(mx, line, op);
      }
    }
    return;
  }

  if( op<0 )  // reverse then repeat
  {
    for(int z=0; z>op; --z)
    {
      for( mx_line_members::reverse_iterator it = members.rbegin()
         ; it!=members.rend(); ++it )
      {
        mx_line_member member = it->first;
        int op = it->second;
        member.interpret(mx, line, op);
      }
    }
    return;
  }
}

void mx_line_member::interpret(MadX const & mx, MadX_line & line, int op)
{
  if( op==0 ) return;

  if( tag==MX_LINE_MEMBER_NAME )
  {
    // member is a name/reference
    string name = any_cast<string>(member);
    MadX_entry_type type = mx.entry_type(name);

    // does the name refer to a madx command?
    if( type==ENTRY_COMMAND )
    { 
      // ok it is a command
      if( mx.command(name, true).is_element() )
      { 
        // now it is a real element
        if( op!=1 )
          throw runtime_error("Line op only applies on sublines, not on elements!");

        // push to the line
        line.insert_element( name );
      }
      else
      {
        throw runtime_error("Line member '" + name + "' is not an element");
      }
    } 
    // or the name referes to a pre-exisitng line
    else if( type==ENTRY_LINE )
    {
      MadX_line const & subline = mx.line(name);
      size_t ne = subline.element_count();
      int repeat = (op>0) ? op : -op;

      for(int z=0; z<repeat; ++z)
      {
        for(size_t i=0; i<ne; ++i)
          line.insert_element( subline.element_name( (op>0) ? (i) : (ne-1-i) ) );
      }
    }
    // something we dont support
    else
    {
      throw runtime_error("Line member '" + name + "' does not exist or not correct type");
    }
  }
  else
  {
    // if it is not a name, it must be a seq
    any_cast<mx_line_seq>(member).interpret(mx, line, op);
  }

  return;
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

void mx_command::interpret(MadX & mx) 
{
  // first execute the command if needed
  execute(mx);

  // push the command into MadX object
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

    // prepare MadX_command object
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

  return;
}

void mx_command::execute(MadX & mx)
{
  if( keyword_ == "call" )
  {
    // pull in the sub-file
    for( attrs_t::const_iterator it = attrs_.begin()
       ; it != attrs_.end(); ++it )
    {
      if( it->name() == "file" )
      {
        string fname = boost::any_cast<string>(it->value());
        mx_tree subroutine;
        parse_int_madx_file(fname, subroutine);
        subroutine.interpret(mx);
        return;
      }
    }
    throw runtime_error("Error executing command 'call'");
  } 
  else if( keyword_ == "sequence" )
  {
    // the "refer" attribute is an unquoted string, not a reference
    for( attrs_t::iterator it = attrs_.begin()
       ; it != attrs_.end(); ++it )
    {
      if( it->name() == "refer" )
      {
        if( it->value().type() == typeid(string) )
        {
          // do nothing
        }
        else if( it->value().type() == typeid(mx_expr) )
        {
          try 
          {
            mx_expr ex = any_cast<mx_expr>(it->value());
            ex = get<nop_t>(get<nop_t>(get<nop_t>(ex).expr).expr).expr;
            string val = boost::get<string>( ex );
            it->set_attr("refer", val);
          } 
          catch(...) 
          {
            throw runtime_error("The 'refer' attribute of sequence '" 
                + label_ + "' is not a string");
          }
        }
        else
        {
          throw runtime_error("The 'refer' attribute of sequence '" 
              + label_ + "' is not a string");
        }
      }

      if( it->name() == "refpos" )
      {
        if( it->value().type() == typeid(string) )
        {
          // do nothing
        }
        else if( it->value().type() == typeid(mx_expr) )
        {
          try 
          {
            mx_expr ex = any_cast<mx_expr>(it->value());
            ex = get<nop_t>(get<nop_t>(get<nop_t>(ex).expr).expr).expr;
            string val = boost::get<string>( ex );
            it->set_attr("refpos", val);
          } 
          catch(...) 
          {
            throw runtime_error("The 'refpos' attribute of sequence '" 
                + label_ + "' is not a string");
          }
        }
        else
        {
          throw runtime_error("The 'refpos' attribute of sequence '" 
              + label_ + "' is not a string");
        }
 
      }
    }
  }
  else if( keyword_ == "beam" )
  {
    // beam can always be referenced with 'beam'
    if( !labeled_ )
    {
      label_ = "beam";
      labeled_ = true;
    }

    // build attributes of beam
    double mass = 0, charge = 0, energy = 0, pc = 0, gamma = 0;
    for ( attrs_t::const_iterator it = attrs_.begin()
        ; it != attrs_.end(); ++it ) 
    {
      if ( it->name() == "particle" ) 
      {
        string particle = any_cast<mx_keyword>( it->value() ).name;
        if (particle == "proton") {
          mass = pconstants::mp;
          charge = pconstants::proton_charge;
        } else if (particle == "antiproton") {
          mass = pconstants::mp;
          charge = pconstants::antiproton_charge;
        } else if (particle == "electron") {
          mass = pconstants::me;
          charge = pconstants::electron_charge;
        } else if (particle == "positron") {
          mass = pconstants::me;
          charge = pconstants::positron_charge;
        } else if (particle == "negmuon") {
          mass = pconstants::mmu;
          charge = pconstants::muon_charge;
        } else if (particle == "posmuon") {
          mass = pconstants::mmu;
          charge = pconstants::antimuon_charge;
        } else {
          throw runtime_error("Unknown particle type " + particle);
        }
      } 
      else if ( it->name() == "mass") 
      {
        mx_expr e = boost::any_cast<mx_expr>( it->value() );
        mass = boost::apply_visitor(mx_calculator(mx), e);
      } 
      else if ( it->name() == "charge") 
      {
        mx_expr e = boost::any_cast<mx_expr>( it->value() );
        charge = boost::apply_visitor(mx_calculator(mx), e);
      } 
      else if ( it->name() == "energy") 
      {
        mx_expr e = boost::any_cast<mx_expr>( it->value() );
        energy = boost::apply_visitor(mx_calculator(mx), e);
      } 
      else if ( it->name() == "pc") 
      {
        mx_expr e = boost::any_cast<mx_expr>( it->value() );
        pc = boost::apply_visitor(mx_calculator(mx), e);
      } 
      else if ( it->name() == "gamma") 
      {
        mx_expr e = boost::any_cast<mx_expr>( it->value() );
        gamma = boost::apply_visitor(mx_calculator(mx), e);
      }
    }

    // makes no change if the particle type is absent
    if( mass==0 ) return;

    Four_momentum four_momentum(mass);

    if (energy > 0) four_momentum.set_total_energy(energy);
    if (pc > 0)     four_momentum.set_momentum(pc);
    if (gamma > 0)  four_momentum.set_gamma(gamma);

    mx_attr attr;
    if (energy == 0) 
    {
      attr.set_attr( "energy", mx_expr(four_momentum.get_total_energy()) );
      ins_attr(attr);
    }
    if (pc == 0) 
    {
      attr.set_attr( "pc", mx_expr(four_momentum.get_momentum()) );
      ins_attr(attr);
    }
    if (gamma == 0) 
    {
      attr.set_attr( "gamma", mx_expr(four_momentum.get_gamma()) );
      ins_attr(attr);
    }
  }
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
bool mx_if_block::evaluate_logic(MadX const & mx) const
{
  return logic_expr.evaluate(mx);
}

void mx_if_block::interpret_block(MadX & mx)
{
  block.interpret(mx);
}

void mx_if_block::print_logic() const
{
  //cout << logic_expr;
}

void mx_if_block::print_block() const
{
  block.print();
}

// if-elseif-else
void mx_if::assign_if(mx_logic const & logic, mx_tree const & block)
{
  if_ = mx_if_block(logic, block);
}

void mx_if::assign_elseif(mx_logic const & logic, mx_tree const & block)
{
  elseif_.push_back(mx_if_block(logic, block));
}

void mx_if::assign_else(mx_tree const & block)
{
  else_ = mx_if_block(true, block);
}

void mx_if::interpret(MadX & mx)
{
  if( !if_.valid() )
    throw runtime_error("mx_if::interpret() Invalid if command");

  if( if_.evaluate_logic(mx) )
  {
    if_.interpret_block(mx);
  }
  else if( !elseif_.empty() )
  {
    for( mx_if_block_v::iterator it = elseif_.begin();
         it != elseif_.end();
         ++it )
    {
      if( it->evaluate_logic(mx) )
      {
        it->interpret_block(mx);
      }
    }
  }
  else if( else_.valid() )
  {
    else_.interpret_block(mx);
  }

  return;
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
void mx_while::assign(mx_logic const & logic, mx_tree const & block)
{
  while_ = mx_if_block(logic, block);
}

void mx_while::interpret(MadX & mx)
{
  return;
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

void mx_statement::assign(mx_line const & st)
{
  value = any(st);
  type = MX_LINE;
}

void mx_statement::interpret(MadX & mx) 
{
  switch( type )
  {
  case MX_COMMAND:
    boost::any_cast<mx_command>(value).interpret(mx);
    break;
  case MX_LINE:
    boost::any_cast<mx_line>(value).interpret(mx);
    break;
  case MX_IF:
    boost::any_cast<mx_if>(value).interpret(mx);
    break;
  case MX_WHILE:
    boost::any_cast<mx_while>(value).interpret(mx);
    break;
  case MX_NULL:
    break;
  default:
    throw runtime_error("mx_statement::interpret()  Unknown statement type");
  }
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
void mx_tree::interpret(MadX & mx) 
{
  for( mx_statements_t::iterator it = statements.begin()
     ; it != statements.end(); ++it )
  {
    it->interpret(mx);
  }
}

void mx_tree::print() const
{
  for_each(statements.begin(), statements.end(), detail::print<mx_statement>);
}


















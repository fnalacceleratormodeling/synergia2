
#include "madx.h"
#include "mx_expr.h"

#include <cmath>
#include <stdexcept>
#include <iostream>


using namespace synergia;
using namespace std;


//===========================================================================
// Helper functions

namespace
{
  string_t
    retrieve_string_from_map( value_map_t const & m
                            , string_t const & k )
  {
    string_t key(k);
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);

    value_map_t::const_iterator it = m.find(key);
    if( it!=m.end() )
    {
      return boost::any_cast<string_t>(it->second.value);
    }
    else
    {
      throw std::runtime_error( "cannot find attribute with name " + key);
    }
  }

  double
    retrieve_number_from_map( value_map_t const & m
                            , string_t const & k
                            , MadX const & global )
  {
    string_t key(k);
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);

    value_map_t::const_iterator it = m.find(key);
    if( it!=m.end() )
    {
      if( it->second.type == NUMBER )
      {
        mx_expr e = boost::any_cast<mx_expr>(it->second.value);
        return boost::apply_visitor(mx_calculator(global), e);
      }
      else
      {
        throw std::runtime_error( "the requested key '" + k + "' cannot be retrieved as a number" );
      }
    }
    else
    {
      throw std::runtime_error( "cannot find attribute with name " + key);
    }
  }

  std::vector<double>
    retrieve_number_seq_from_map( value_map_t const & m
                                , string_t const & k
                                , MadX const & global )
  {
    string_t key(k);
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);

    value_map_t::const_iterator it = m.find(key);
    if( it!=m.end() )
    {
     if( it->second.type == ARRAY )
      {
        mx_exprs es = boost::any_cast<mx_exprs>(it->second.value);
        std::vector<double> vd;

        for( mx_exprs::const_iterator it = es.begin()
           ; it != es.end(); ++it )
        {
          vd.push_back( boost::apply_visitor(mx_calculator(global), *it) );
        }
        return vd;
      }
      else
      {
        throw std::runtime_error( "the requested key '" + k + "' cannot be retrieved as a sequence of number" );
      }
    }
    else
    {
      throw std::runtime_error( "cannot find attribute with name " + key);
    }
  }

  MadX_command
    resolve_command(MadX_command const & cmd, MadX const & mx, bool resolve)
  {
    if( !resolve || !cmd.is_reference() )  return cmd;

    MadX_command result = mx.command(cmd.name());
    result.merge_with_overwrite(cmd);
    return result;
  }

}



//===========================================================================
// MadX_command

string_t
  MadX_command::name() const
{
  return name_;
}

string_t
  MadX_command::label() const
{
  return label_;
}

size_t
  MadX_command::attribute_count() const
{
  return attributes_.size();
}

std::vector<string_t>
  MadX_command::attribute_names() const
{
  std::vector<string_t> names;
  for( value_map_t::const_iterator it = attributes_.begin()
     ; it!= attributes_.end(); ++it )
    names.push_back(it->first);
  return names;
}

MadX_value_type
  MadX_command::attribute_type( string_t const & name ) const
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  value_map_t::const_iterator it = attributes_.find(key);
  if( it!=attributes_.end() )
  {
    return it->second.type;
  }
  else
  {
    throw std::runtime_error( "MadX_command::attribute_type:"
        " cannot find attribute with name " + key);
  }
}

string_t
  MadX_command::attribute_as_string( string_t const & name ) const
{
  return retrieve_string_from_map(attributes_, name);
}

double
  MadX_command::attribute_as_number( string_t const & name ) const
{
  return retrieve_number_from_map(attributes_, name, *mx);
}

bool
  MadX_command::attribute_as_boolean( string_t const & name ) const
{
  return std::abs(retrieve_number_from_map(attributes_, name, *mx)) > 1e-10;
}

std::vector<double>
  MadX_command::attribute_as_number_seq( string_t const & name ) const
{
  return retrieve_number_seq_from_map(attributes_, name, *mx);
}

void
  MadX_command::set_parent( MadX const & parent )
{
  mx = &parent;
}

void
  MadX_command::set_name( string_t const & name, MadX_command_type type )
{
  name_ = name;
  std::transform(name_.begin(), name_.end(), name_.begin(), ::tolower);
  type_ = type;
}

void
  MadX_command::set_label( string_t const & label)
{
  label_ = label;
}

MadX_command_type
  MadX_command::type() const
{
  return type_;
}

bool
  MadX_command::is_element() const
{
  return type_ == ELEMENT;
}

bool
  MadX_command::is_reference() const
{
  return type_ == ELEMENT_REF;
}

bool
  MadX_command::is_command() const
{
  return type_ == EXECUTABLE;
}

void
  MadX_command::insert_attribute( string_t const & name, string_t const & value )
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  MadX_value v;
  v.value = boost::any(value);
  v.type  = value.empty() ? NONE : STRING;

  attributes_.insert(std::make_pair(key, v));
}

void
  MadX_command::insert_attribute( string_t const & name, mx_expr const & value )
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  MadX_value v;
  v.value = boost::any(value);
  v.type  = NUMBER;

  attributes_.insert(std::make_pair(key, v));
}

void
  MadX_command::insert_attribute( string_t const & name, mx_exprs const & value )
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  MadX_value v;
  v.value = boost::any(value);
  v.type  = ARRAY;

  attributes_.insert(std::make_pair(key, v));
}

void
  MadX_command::merge_with_overwrite(MadX_command const & other)
{
  value_map_t map(other.attributes_);
  map.insert(attributes_.begin(), attributes_.end());
  std::swap(attributes_, map);
}

void
  MadX_command::merge(MadX_command const & other)
{
  attributes_.insert(other.attributes_.begin(), other.attributes_.end());
}


//===========================================================================
// MadX_line
size_t
  MadX_line::element_count() const
{
  return elements_.size();
}

string_t
  MadX_line::element_name(size_t idx) const
{
  return elements_[idx];
}

MadX_command
  MadX_line::element(size_t idx, bool resolve) const
{
  return parent.command(elements_[idx], resolve);
}

void
  MadX_line::insert_element(string_t const & ele)
{
  string_t e(ele);
  std::transform(e.begin(), e.end(), e.begin(), ::tolower);

  elements_.push_back(e);
}

void
  MadX_line::print() const
{
  cout << "( ";
  for( vector<string>::const_iterator it=elements_.begin()
     ; it!=elements_.end(); ++it )
  {
    cout << *it;
    if( it+1 != elements_.end() ) cout << ", ";
  }
  cout << " )\n";
}

//===========================================================================
// MadX_sequence

string_t
  MadX_sequence::label() const
{
  return lbl;
}

double
  MadX_sequence::length() const
{
  return l;
}

size_t
  MadX_sequence::element_count() const
{
  return seq_.size();
}

MadX_command
  MadX_sequence::element(size_t idx, bool resolve) const
{
  return resolve_command(seq_[idx], parent, resolve);
}

void
  MadX_sequence::set_label(string_t const & label)
{
  lbl = label;
}

void
  MadX_sequence::set_length(double length)
{
  l = length;
}

void
  MadX_sequence::add_element(MadX_command const & cmd)
{
  seq_.push_back(cmd);
}

void
  MadX_sequence::reset()
{
  lbl = string_t();
  l = 0.0;
  seq_.clear();
}

void
  MadX_sequence::print() const
{

}


//===========================================================================
// MadX

string_t
  MadX::variable_as_string( string_t const & name ) const
{
  return retrieve_string_from_map( variables_, name );
}

double
  MadX::variable_as_number( string_t const & name ) const
{
  return retrieve_number_from_map( variables_, name, *this );
}

bool
  MadX::variable_as_boolean( string_t const & name ) const
{
  return std::abs(retrieve_number_from_map( variables_, name, *this )) > 1e-10;
}

std::vector<double>
  MadX::variable_as_number_seq( string_t const & name ) const
{
  return retrieve_number_seq_from_map( variables_, name, *this );
}

size_t
  MadX::command_count() const
{
  return cmd_seq_.size();
}

std::vector<string_t >
  MadX::commands() const
{
    std::vector<string_t > commands;
    for(commands_v_t::const_iterator it = cmd_seq_.begin();
            it != cmd_seq_.end(); ++it)
    {
        commands.push_back(it->name());
    }
    return commands;
}

MadX_command
  MadX::command( size_t idx, bool resolve ) const
{
  return resolve_command(cmd_seq_[idx], *this, resolve);
}

size_t
  MadX::label_count() const
{
  return cmd_map_.size();
}

std::vector<string_t>
  MadX::command_labels() const
{
  std::vector<string_t> labels;
  for( commands_m_t::const_iterator it = cmd_map_.begin()
     ; it!=cmd_map_.end(); ++it )
    labels.push_back(it->first);
  return labels;
}

MadX_command
  MadX::command( string_t const & label, bool resolve ) const
{
  string_t key(label);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  commands_m_t::const_iterator it = cmd_map_.find(key);

  if( it!=cmd_map_.end() )
  {
    return resolve_command(it->second, *this, resolve);
  }
  else
  {
    throw std::runtime_error( "cannot find command with label " + key );
  }
}

size_t
  MadX::line_count() const
{
  return lines_.size();
}

std::vector<string_t>
  MadX::line_labels() const
{
  std::vector<string_t> labels;
  for( lines_m_t::const_iterator it = lines_.begin()
     ; it!=lines_.end(); ++it)
    labels.push_back(it->first);
  return labels;
}

MadX_line const &
  MadX::line( string_t const & line ) const
{
  string_t key(line);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  lines_m_t::const_iterator it = lines_.find(key);

  if( it!=lines_.end() )
    return it->second;
  else
    throw std::runtime_error( "cannot find line with label " + line );
}

size_t
  MadX::sequence_count() const
{
  return seqs_.size();
}

std::vector<string_t>
  MadX::sequence_labels() const
{
  std::vector<string_t> labels;
  for( sequences_m_t::const_iterator it = seqs_.begin()
     ; it!=seqs_.end(); ++it)
    labels.push_back(it->first);
  return labels;
}

MadX_sequence const &
  MadX::sequence( string_t const & seq ) const
{
  string_t key(seq);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  sequences_m_t::const_iterator it = seqs_.find(key);

  if( it!=seqs_.end() )
    return it->second;
  else
    throw std::runtime_error( "cannot find sequence with label " + seq );
}

MadX_sequence const &
  MadX::current_sequence( ) const
{
  return cur_seq_;
}

MadX_sequence &
  MadX::current_sequence( )
{
  return cur_seq_;
}

MadX_entry_type
  MadX::entry_type(string_t const & entry) const
{
  string_t key(entry);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  if( variables_.find(key) != variables_.end() )
    return ENTRY_VARIABLE;

  if( cmd_map_.find(key) != cmd_map_.end() )
    return ENTRY_COMMAND;

  if( lines_.find(key) != lines_.end() )
    return ENTRY_LINE;

  if( seqs_.find(key) != seqs_.end() )
    return ENTRY_SEQUENCE;

  return ENTRY_NULL;
}

void
  MadX::insert_variable(string_t const & name, string_t const & value)
{
  string_t key(name);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  MadX_value v;
  v.value = boost::any(value);
  v.type  = value.empty() ? NONE : STRING;

  variables_.insert(std::make_pair(key, v));
}

void
  MadX::insert_variable(string_t const & name, mx_expr const & value)
{
  string_t key(name);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  MadX_value v;
  v.value = boost::any(value);
  v.type  = NUMBER;

  variables_.insert( std::make_pair(key, v) );
}

void
  MadX::insert_variable(string_t const & name, mx_exprs const & value)
{
  string_t key(name);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  MadX_value v;
  v.value = boost::any(value);
  v.type  = ARRAY;

  variables_.insert( std::make_pair(key, v) );
}

void
  MadX::insert_label(string_t const & name, MadX_command const & cmd)
{
  string_t key(name);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  cmd_map_.insert( std::make_pair(key, cmd) );
  cmd_map_[key].set_parent(*this);

  // TODO : execution part should be moved to mx_tree
  execute_command(key, cmd);
}

void
  MadX::insert_line(string_t const & name, MadX_line const & line)
{
  string_t key(name);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  lines_.insert( std::make_pair(key, line) );
}

void
  MadX::insert_command(MadX_command const & cmd)
{
  cmd_seq_.push_back(cmd);
  cmd_seq_.back().set_parent(*this);

  // TODO : execution part should be moved to mx_tree
  execute_command("", cmd);
}

void
  MadX::start_sequence()
{
  building_seq_ = true;
}

void
  MadX::end_sequence()
{
  building_seq_ = false;
}

void
  MadX::execute_command(string_t const & label, MadX_command const & cmd)
{
  if( cmd.name()=="sequence" )
  {
    // TODO: start building sequence
    building_seq_ = true;
    cur_seq_.set_label( label );
    cur_seq_.set_length( cmd.attribute_as_number("l") );
  }
  else if( cmd.name()=="endsequence" )
  {
    // TODO: finish building the sequence and push it to the madx object
    building_seq_ = false;
    seqs_.insert(std::make_pair(cur_seq_.label(), cur_seq_));
    cur_seq_.reset();
  }
  else if( building_seq_ )
  {
    // TODO: command must have "at" field
    // throw if no 'at'

    // push to cur_seq_ object
    cur_seq_.add_element( cmd );
  }
}

void
  MadX::print() const
{
  value_map_t::const_iterator it = variables_.begin();
  for(; it!=variables_.end(); ++it )
  {
    std::string key = it->first;
    double val = variable_as_number(key);
    std::cout << key << " = " << val << "\n";
  }

  for( lines_m_t::const_iterator it=lines_.begin()
     ; it!=lines_.end(); ++it )
  {
    cout << it->first << " : line = ";
    it->second.print();
  }
}





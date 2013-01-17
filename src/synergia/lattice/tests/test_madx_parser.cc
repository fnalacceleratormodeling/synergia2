#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

#include <iostream>

#include "synergia/lattice/mx_parse.h"

using namespace std;
using namespace synergia;

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct)
{
}

BOOST_AUTO_TEST_CASE(blank_line)
{
  string str = "";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
}

BOOST_AUTO_TEST_CASE(garbage)
{
  string str = "your momma uses portions->of madx syntax";
  MadX   mx;

  BOOST_CHECK( !parse_madx( str, mx ) );
}

#if 0
BOOST_AUTO_TEST_CASE(comment1)
{
  string str = "//your momma uses portions->of madx syntax";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
}

BOOST_AUTO_TEST_CASE(comment2)
{
  string str = "!your momma uses portions->of madx syntax";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
}
#endif

BOOST_AUTO_TEST_CASE(variable_assignment)
{
  string str = "x=1;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.variable_as_number("x"), 1 );
}


BOOST_AUTO_TEST_CASE(mod_variable_assignment)
{
  string str = "x:=1;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.variable_as_number("x"), 1 );
}

BOOST_AUTO_TEST_CASE(newline_separation)
{
  string str = "x=1;\n y=2;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.variable_as_number("x"), 1 );
  BOOST_CHECK_EQUAL( mx.variable_as_number("y"), 2 );
}

BOOST_AUTO_TEST_CASE(semicolon_separation)
{
  string str = "x=1;y=2;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.variable_as_number("x"), 1 );
  BOOST_CHECK_EQUAL( mx.variable_as_number("y"), 2 );
}

BOOST_AUTO_TEST_CASE(variable_assignment_expression)
{
  string str = "foo.bar=pi*sin(1.2e-4)^0.69;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_CLOSE( mx.variable_as_number("foo.bar"),
                     boost::math::constants::pi<double>() * pow(sin(1.2e-4), 0.69),
                     tolerance );
}

BOOST_AUTO_TEST_CASE(caps_variable_assignment)
{
  string str = "X=1;Y=X;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.variable_as_number("y"), 1 );
}

// ==========================================================================
// commands

BOOST_AUTO_TEST_CASE(command)
{
  string str = "beam;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.command_count(), 1 );

  MadX_command cmd = mx.command(0);
  BOOST_CHECK_EQUAL( cmd.name(), "beam" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 0 );
}

#if 0
BOOST_AUTO_TEST_CASE(upper_command)
{
  string str = "BEAM;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.command_count(), 1 );

  MadX_command cmd = mx.command(0);
  BOOST_CHECK_EQUAL( cmd.name(), "beam" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 0 );
}
#endif

BOOST_AUTO_TEST_CASE(command_attrs)
{
  string str = "beam,a=1,B=3*(4+5);";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.command_count(), 1 );

  MadX_command cmd = mx.command(0);
  BOOST_CHECK_EQUAL( cmd.name(), "beam" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 2 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_number("a"), 1 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_number("b"), 3*(4+5) );
}

BOOST_AUTO_TEST_CASE(command_str_attrs1)
{
  string str = "title, S = \"Tevatron Collider Run II Lattice\";";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.command_count(), 1 );

  MadX_command cmd = mx.command(0);
  BOOST_CHECK_EQUAL( cmd.name(), "title" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 1 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_string("s"), "Tevatron Collider Run II Lattice" );
}

#if 0
BOOST_AUTO_TEST_CASE(command_str_attrs2)
{
  string str = "title, S = 'Tevatron Collider Run II Lattice';";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.command_count(), 1 );

  MadX_command cmd = mx.command(0);
  BOOST_CHECK_EQUAL( cmd.name(), "title" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 1 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_string("s"), "Tevatron Collider Run II Lattice" );
}
#endif

#if 0
BOOST_AUTO_TEST_CASE(command_particle_attrs)
{
  string str = "beam, particle=proton;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.command_count(), 1 );

  MadX_command cmd = mx.command(0);
  BOOST_CHECK_EQUAL( cmd.name(), "beam" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 1 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_string("particle"), "proton");
}
#endif

BOOST_AUTO_TEST_CASE(command_assign)
{
  string str = "q1: quadrupole,l=3.14;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("q1");
  BOOST_CHECK_EQUAL( cmd.name(), "quadrupole" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 1 );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("l"), 3.14, tolerance );
}

BOOST_AUTO_TEST_CASE(subscripted_ident)
{
  string str = "foo: quadrupole, a=1; x=foo->a;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.variable_as_number("x"), 1 );
}

BOOST_AUTO_TEST_CASE(continuation)
{
  string str = "q1: quadrupole,l=\n3.14,k1=0.2;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("q1");
  BOOST_CHECK_EQUAL( cmd.name(), "quadrupole" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 2 );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("l"), 3.14, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("k1"), 0.2, tolerance );
}

BOOST_AUTO_TEST_CASE(continuation2)
{
  string str = "q2: quadrupole,l=! )junk)\n3.14,k1=0.2;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("q2");
  BOOST_CHECK_EQUAL( cmd.name(), "quadrupole" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 2 );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("l"), 3.14, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("k1"), 0.2, tolerance );
}

BOOST_AUTO_TEST_CASE(continuation3)
{
  string str = "q3: quadrupole,l=3.14,\nk1=0.2;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("q3");
  BOOST_CHECK_EQUAL( cmd.name(), "quadrupole" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 2 );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("l"), 3.14, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("k1"), 0.2, tolerance );
}

#if 0
BOOST_AUTO_TEST_CASE(continuation4)
{
  string str = "twiss, save,   betx=28.871,alfx=-0.069,mux=0.0,dx=2.682,dpx=-0.073,\nbety= 5.264,alfy=-0.006,muy=0.0,dy=0.0,dpy=0.0;";
  MadX   mx;

  BOOST_CHECK( parse_madx( str, mx ) );
}
#endif


#if 0
#include <boost/any.hpp>

using namespace boost;

int main()
{
  string fname = "x:=sp+3; \n sp=4;";
  synergia::MadX mx;

  //cout << "Please input the lattice file name: ";
  //cin  >> fname;

  //bool result = parse_madx( fname, mx );
  bool result = parse_madx_file("psb_lattice.madx", mx);

  if( result ) cout << "Parsing succeeded!\n";
  else         cout << "Parsing failed!\n";

  //mx.print();

  cout << mx.command("br.qfo12").attribute_as_number("k1") << "\n";
  cout << mx.command("qfo").attribute_as_number("k1") << "\n";

#if 0
  string val = "3+4*12+x";
  any r;
  result = evaluate_madx_value( val, r, &mx );
  cout << any_cast<double>(r) << "\n";

  string mx_str = " a:=1; b=a; if(aaa) {hello;a; } while(bbb){heyday;world;}!nonsense\nx=3;";
  synergia::mx_tree tree;
  result = parse_madx_tree( mx_str, tree );
  tree.print();
#endif
}
#endif

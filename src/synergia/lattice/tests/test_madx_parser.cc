#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

#include <iostream>

#include "synergia/foundation/physical_constants.h"
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

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
}

BOOST_AUTO_TEST_CASE(garbage)
{
  string str = "your momma uses portions->of madx syntax";
  MadX   mx;

  BOOST_CHECK_THROW( parse_madx( str, mx ), std::runtime_error );
}

#if 0
BOOST_AUTO_TEST_CASE(comment1)
{
  string str = "//your momma uses portions->of madx syntax";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
}

BOOST_AUTO_TEST_CASE(comment2)
{
  string str = "!your momma uses portions->of madx syntax";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
}
#endif

BOOST_AUTO_TEST_CASE(variable_assignment)
{
  string str = "x=1;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.variable_as_number("x"), 1 );
}


BOOST_AUTO_TEST_CASE(mod_variable_assignment)
{
  string str = "x:=1;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.variable_as_number("x"), 1 );
}

BOOST_AUTO_TEST_CASE(keyword_leading_values)
{
  // n.b. "pine" starts with "pi"
  string str = "xpine = 3; v = xpine; endmark: marker;";
  str += "fodo: sequence, refer=entry, l=12.0; endmark, at=5.0; endsequence;";
  MadX mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );

  MadX_sequence seq = mx.sequence("fodo");
  BOOST_CHECK_EQUAL( seq.element_count(), 1 );
  BOOST_CHECK_CLOSE( seq.element_at(0), 5.0, tolerance);

  MadX_command cmd = seq.element(0);
  BOOST_CHECK_EQUAL( cmd.name(), "marker" );
}

BOOST_AUTO_TEST_CASE(mad_constants)
{
  string str = "a = pi; b = twopi; c = degrad; d = raddeg; ee = e; "
               "f = emass; g = pmass; h = mumass; i = clight; j = qelect; ";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  
  BOOST_CHECK_CLOSE( mx.variable_as_number("a"), boost::math::constants::pi<double>(),         tolerance);
  BOOST_CHECK_CLOSE( mx.variable_as_number("b"), boost::math::constants::two_pi<double>(),     tolerance);
  BOOST_CHECK_CLOSE( mx.variable_as_number("c"), 180.0 / boost::math::constants::pi<double>(), tolerance);
  BOOST_CHECK_CLOSE( mx.variable_as_number("d"), boost::math::constants::pi<double>() / 180.0, tolerance);
  BOOST_CHECK_CLOSE( mx.variable_as_number("ee"), boost::math::constants::e<double>(),         tolerance);

  BOOST_CHECK_CLOSE( mx.variable_as_number("f"), pconstants::me,  tolerance);
  BOOST_CHECK_CLOSE( mx.variable_as_number("g"), pconstants::mp,  tolerance);
  BOOST_CHECK_CLOSE( mx.variable_as_number("h"), pconstants::mmu, tolerance);
  BOOST_CHECK_CLOSE( mx.variable_as_number("i"), pconstants::c,   tolerance);
  BOOST_CHECK_CLOSE( mx.variable_as_number("j"), pconstants::e,   tolerance);
}



BOOST_AUTO_TEST_CASE(newline_separation)
{
  string str = "x=1;\n y=2;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.variable_as_number("x"), 1 );
  BOOST_CHECK_EQUAL( mx.variable_as_number("y"), 2 );
}

BOOST_AUTO_TEST_CASE(semicolon_separation)
{
  string str = "x=1;y=2;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.variable_as_number("x"), 1 );
  BOOST_CHECK_EQUAL( mx.variable_as_number("y"), 2 );
}

BOOST_AUTO_TEST_CASE(floating_point)
{
  // 7/32 has an exact floating point representation
  string str = "x=1.234;y=0.21875; ze2=.21875e2; wep2=.21875E2; wep02=.21875e+02;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_CLOSE( mx.variable_as_number("x"), 1.234, tolerance);
  BOOST_CHECK_EQUAL( mx.variable_as_number("y"), 0.21875 );
  BOOST_CHECK_EQUAL( mx.variable_as_number("ze2"), 21.875);
  BOOST_CHECK_EQUAL( mx.variable_as_number("wep2"), 21.875);
  BOOST_CHECK_EQUAL( mx.variable_as_number("wep02"), 21.875);
}

BOOST_AUTO_TEST_CASE(variable_assignment_expression)
{
  string str = "foo.bar=pi*sin(1.2e-4)^0.69;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_CLOSE( mx.variable_as_number("foo.bar"),
                     boost::math::constants::pi<double>() * pow(sin(1.2e-4), 0.69),
                     tolerance );
}

BOOST_AUTO_TEST_CASE(caps_variable_assignment)
{
  string str = "X=1;Y=X;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.variable_as_number("y"), 1 );
}

// ==========================================================================
// commands

BOOST_AUTO_TEST_CASE(command)
{
  string str = "beam;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("beam");
  BOOST_CHECK_EQUAL( cmd.name(), "beam" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 0 );
}

#if 0
BOOST_AUTO_TEST_CASE(upper_command)
{
  string str = "BEAM;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
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

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("beam");
  BOOST_CHECK_EQUAL( cmd.name(), "beam" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 2 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_number("a"), 1 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_number("b"), 3*(4+5) );
}

#if 0
BOOST_AUTO_TEST_CASE(command_str_attrs1)
{
  string str = "title, S = \"Tevatron Collider Run II Lattice\";";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.command_count(), 1 );

  MadX_command cmd = mx.command(0);
  BOOST_CHECK_EQUAL( cmd.name(), "title" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 1 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_string("s"), "Tevatron Collider Run II Lattice" );
}

BOOST_AUTO_TEST_CASE(command_str_attrs2)
{
  string str = "title, S = 'Tevatron Collider Run II Lattice';";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.command_count(), 1 );

  MadX_command cmd = mx.command(0);
  BOOST_CHECK_EQUAL( cmd.name(), "title" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 1 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_string("s"), "Tevatron Collider Run II Lattice" );
}
#endif

BOOST_AUTO_TEST_CASE(command_particle_attrs)
{
  string str = "beAM, particle=proton;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("beam");
  BOOST_CHECK_EQUAL( cmd.name(), "beam" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 5 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_string("particle"), "proton");
}

BOOST_AUTO_TEST_CASE(command_special_attrs1)
{
  string str = "mp: multipole, knl:={0, 1, 1}, type=octpn;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("mp");
  BOOST_CHECK_EQUAL( cmd.name(), "multipole" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 2 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_string("type"), "octpn");
}

BOOST_AUTO_TEST_CASE(command_special_attrs2)
{
  string str = "mp: multipole, knl:={0, 1, 1}, TYPE=wgl;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("mp");
  BOOST_CHECK_EQUAL( cmd.name(), "multipole" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 2 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_string("type"), "wgl");
}

BOOST_AUTO_TEST_CASE(command_special_attrs3)
{
  string str = "mp: multipole, knl:={0, 1, 1}, type=\"special\";";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("mp");
  BOOST_CHECK_EQUAL( cmd.name(), "multipole" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 2 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_string("type"), "special");
}

BOOST_AUTO_TEST_CASE(command_omitted_comma)
{
  string str = "call file = './foo.dbx';";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.command_count(), 0 );

  BOOST_CHECK_EQUAL( mx.variable_as_number("a"), 3 );
  BOOST_CHECK_EQUAL( mx.variable_as_number("b"), 2 );
}

BOOST_AUTO_TEST_CASE(command_beam_particle)
{
  string str = " BEAM, PARTICLE=Proton, MASS=0.93827, CHARGE=1., ENERGY=0.93827 + 0.160;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );

  MadX_command cmd = mx.command("beam");
  BOOST_CHECK_EQUAL( cmd.name(), "beam" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 6 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_number("mass"), 0.93827 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_number("charge"), 1 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_number("energy"), 0.93827 + 0.160 );
}

BOOST_AUTO_TEST_CASE(command_beam_particle_abbreviate)
{
  string str = " BEAM, PARTICLE=Prot, MASS=0.93827, CHARGE=1., ENERGY=0.93827 + 0.160;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );

  MadX_command cmd = mx.command("beam");
  BOOST_CHECK_EQUAL( cmd.name(), "beam" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 6 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_number("mass"), 0.93827 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_number("charge"), 1 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_number("energy"), 0.93827 + 0.160 );
}

BOOST_AUTO_TEST_CASE(command_assign)
{
  string str = "q1: quadrupole,l=3.14;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
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

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.variable_as_number("x"), 1 );
}

BOOST_AUTO_TEST_CASE(vector_attribute)
{
    string str = "mpfi1: multipole,knl:={ 0, 1.1, 2.2, 3.3, 4.4, 5.5 };";
    MadX mx;
    BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
    MadX_command cmd = mx.command("mpfi1");
    std::vector<double > knl(cmd.attribute_as_number_seq("knl"));

    BOOST_CHECK_EQUAL(knl.size(), 6);
    const double tolerance = 1.0e-12;
    BOOST_CHECK_CLOSE(knl.at(0), 0, tolerance);
    BOOST_CHECK_CLOSE(knl.at(1), 1.1, tolerance);
    BOOST_CHECK_CLOSE(knl.at(2), 2.2, tolerance);
    BOOST_CHECK_CLOSE(knl.at(3), 3.3, tolerance);
    BOOST_CHECK_CLOSE(knl.at(4), 4.4, tolerance);
    BOOST_CHECK_CLOSE(knl.at(5), 5.5, tolerance);
}

BOOST_AUTO_TEST_CASE(continuation)
{
  string str = "q1: quadrupole,l=\n3.14,k1=0.2;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
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

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
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

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("q3");
  BOOST_CHECK_EQUAL( cmd.name(), "quadrupole" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 2 );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("l"), 3.14, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("k1"), 0.2, tolerance );
}

BOOST_AUTO_TEST_CASE(tkicker)
{
  string str = "tk: tkicker, hkick=0.01;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("tk");
  BOOST_CHECK_EQUAL( cmd.name(), "tkicker" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 1 );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("hkick"), 0.01, tolerance );
}


BOOST_AUTO_TEST_CASE(line_expansion)
{
  string str = "a:rbend; b:rbend; c:rbend; d:rbend; e:rbend; f:rbend; g:rbend; h:rbend;"
	"r:line=(g,h); s:line=(c,r,d); t:line=(2*s,2*(e,f),-s,-(a,b));";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.line_count(), 3 );

  MadX_line line = mx.line("t");
  BOOST_CHECK_EQUAL( line.element_count(), 18 );
  BOOST_CHECK_EQUAL( line.element_name(0), "c" );
  BOOST_CHECK_EQUAL( line.element_name(1), "g" );
  BOOST_CHECK_EQUAL( line.element_name(2), "h" );
  BOOST_CHECK_EQUAL( line.element_name(3), "d" );
  BOOST_CHECK_EQUAL( line.element_name(4), "c" );
  BOOST_CHECK_EQUAL( line.element_name(5), "g" );
  BOOST_CHECK_EQUAL( line.element_name(6), "h" );
  BOOST_CHECK_EQUAL( line.element_name(7), "d" );
  BOOST_CHECK_EQUAL( line.element_name(8), "e" );
  BOOST_CHECK_EQUAL( line.element_name(9), "f" );
  BOOST_CHECK_EQUAL( line.element_name(10), "e" );
  BOOST_CHECK_EQUAL( line.element_name(11), "f" );
  BOOST_CHECK_EQUAL( line.element_name(12), "d" );
  BOOST_CHECK_EQUAL( line.element_name(13), "h" );
  BOOST_CHECK_EQUAL( line.element_name(14), "g" );
  BOOST_CHECK_EQUAL( line.element_name(15), "c" );
  BOOST_CHECK_EQUAL( line.element_name(16), "b" );
  BOOST_CHECK_EQUAL( line.element_name(17), "a" );
}

BOOST_AUTO_TEST_CASE(matrix)
{
  string str = "m1: matrix, type=abc, L=1.2, kick1=0.001, kick2=0.001, kick3=0.002, kick4=0.002, kick5=0.003, kick6=0.003, rm11=1.1, rm12=1.2, rm32=3.2, rm54=5.4, tm111=1.11, tm321=3.21;";

  MadX mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
  BOOST_CHECK_EQUAL( mx.label_count(), 1 );

  MadX_command cmd = mx.command("m1");
  BOOST_CHECK_EQUAL( cmd.name(), "matrix" );
  BOOST_CHECK_EQUAL( cmd.attribute_count(), 14 );
  BOOST_CHECK_EQUAL( cmd.attribute_as_string("type"), "abc" );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("l"), 1.2, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("kick1"), 0.001, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("kick2"), 0.001, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("kick3"), 0.002, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("kick4"), 0.002, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("kick5"), 0.003, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("kick6"), 0.003, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("rm11"), 1.1, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("rm12"), 1.2, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("rm32"), 3.2, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("rm54"), 5.4, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("tm111"), 1.11, tolerance );
  BOOST_CHECK_CLOSE( cmd.attribute_as_number("tm321"), 3.21, tolerance );
}

#if 0
BOOST_AUTO_TEST_CASE(continuation4)
{
  string str = "twiss, save,   betx=28.871,alfx=-0.069,mux=0.0,dx=2.682,dpx=-0.073,\nbety= 5.264,alfy=-0.006,muy=0.0,dy=0.0,dpy=0.0;";
  MadX   mx;

  BOOST_CHECK_NO_THROW( parse_madx( str, mx ) );
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

#include "synergia/utils/catch.hpp"

#include <iostream>
#include <boost/math/constants/constants.hpp>

#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/mx_parse.h"

using namespace std;
using namespace synergia;

const double tolerance = 1.0e-12;

TEST_CASE("line_expansion")
{
    using namespace synergia;

    std::string str = R"(
        a: rbend; 
        b: rbend; 
        c: rbend; 
        d: rbend; 
        e: rbend; 
        f: rbend; 
        g: rbend; 
        h: rbend;
	    r: line=(g,h); 
        s: line=(c,r,d); 
        t: line=(,,, 2*s 2*(e,f),,-s,,,, ,-(a,b),,);
    )";

    MadX   mx;

    REQUIRE_NOTHROW( parse_madx( str, mx ) );
    REQUIRE( mx.line_count() == 3 );

    MadX_line line = mx.line("t");

    CHECK( line.element_count() == 18 );
    CHECK( line.element_name(0) == "c" );
    CHECK( line.element_name(1) == "g" );
    CHECK( line.element_name(2) == "h" );
    CHECK( line.element_name(3) == "d" );
    CHECK( line.element_name(4) == "c" );
    CHECK( line.element_name(5) == "g" );
    CHECK( line.element_name(6) == "h" );
    CHECK( line.element_name(7) == "d" );
    CHECK( line.element_name(8) == "e" );
    CHECK( line.element_name(9) == "f" );
    CHECK( line.element_name(10) == "e" );
    CHECK( line.element_name(11) == "f" );
    CHECK( line.element_name(12) == "d" );
    CHECK( line.element_name(13) == "h" );
    CHECK( line.element_name(14) == "g" );
    CHECK( line.element_name(15) == "c" );
    CHECK( line.element_name(16) == "b" );
    CHECK( line.element_name(17) == "a" );
}


TEST_CASE("line_expansion2")
{
    using namespace synergia;

    std::string str = R"(
        bpm: drift, l=0.01;

        fq1: quadrupole, l=1.0, k1=0.02;
        d1: drift, l=1.0;
        dq2: quadrupole, l=1.0, k1=-0.02;
        d2: drift, l=1.0;
        d3: drift, l=1.0;
        bpm2: bpm;

        line1 : line=(fq1, bpm1, d1);
        bpm1: bpm;

        line2: line=(dq2, d2, bpm1, d3);

        trouble: line=(line1, line2);
    )";

    MadX   mx;

    REQUIRE_NOTHROW( parse_madx( str, mx ) );

    MadX_line line = mx.line("trouble");
    //line.print();

    CHECK( line.element_count() == 7 );
    CHECK( line.element_name(0) == "fq1" );
    CHECK( line.element_name(1) == "bpm1" );
    CHECK( line.element_name(2) == "d1" );
    CHECK( line.element_name(3) == "dq2" );
    CHECK( line.element_name(4) == "d2" );
    CHECK( line.element_name(5) == "bpm1" );
    CHECK( line.element_name(6) == "d3" );
}

TEST_CASE("foil element")
{
    using namespace synergia;

    std::string str = R"(
        f: foil, xmin=-0.05, xmax=0.05, ymin=-0.03, ymax=0.03, thick=600;
    )";

    MadX   mx;

    REQUIRE_NOTHROW( parse_madx( str, mx ) );

    auto foil = mx.command("f");

    CHECK(foil.is_element());
    CHECK(!foil.is_reference());
    CHECK(!foil.is_command());

    CHECK(foil.attribute_as_number("xmin") == Approx(-0.05).margin(tolerance));
    CHECK(foil.attribute_as_number("xmax") == Approx( 0.05).margin(tolerance));
    CHECK(foil.attribute_as_number("ymin") == Approx(-0.03).margin(tolerance));
    CHECK(foil.attribute_as_number("ymax") == Approx( 0.03).margin(tolerance));
    CHECK(foil.attribute_as_number("thick") == Approx(600).margin(tolerance));
}


#if 1
TEST_CASE("construct")
{
  MadX mx;
}

TEST_CASE("blank_line")
{
  string str = "";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
}

TEST_CASE("garbage")
{
  string str = "your momma uses portions->of madx syntax";
  MadX   mx;

  CHECK_THROWS_AS( parse_madx( str, mx ), std::runtime_error );
}

TEST_CASE("comment1")
{
  string str = "//your momma uses portions->of madx syntax";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
}

TEST_CASE("comment2")
{
  string str = "!your momma uses portions->of madx syntax";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
}

TEST_CASE("variable_assignment")
{
  string str = "x=1;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.variable_as_number("x") == 1 );
}


TEST_CASE("mod_variable_assignment")
{
  string str = "x:=1;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.variable_as_number("x") == 1 );
}

TEST_CASE("keyword_leading_values")
{
  // n.b. "pine" starts with "pi"
  string str = "xpine = 3; v = xpine; endmark: marker;";
  str += "fodo: sequence, refer=entry, l=12.0; endmark, at=5.0; endsequence;";
  MadX mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );

  MadX_sequence seq = mx.sequence("fodo");
  CHECK( seq.element_count() == 1 );
  CHECK( seq.element_at(0) == Approx(5.0).margin(tolerance));

  MadX_command cmd = seq.element(0);
  CHECK( cmd.name() == "marker" );
}

TEST_CASE("mad_constants")
{
  string str = "a = pi; b = twopi; c = degrad; d = raddeg; ee = e; "
               "f = emass; g = pmass; h = mumass; i = clight; j = qelect; ";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  
  CHECK( mx.variable_as_number("a") == Approx(boost::math::constants::pi<double>()).margin(tolerance));
  CHECK( mx.variable_as_number("b") == Approx(boost::math::constants::two_pi<double>()).margin(tolerance));
  CHECK( mx.variable_as_number("c") == Approx(180.0 / boost::math::constants::pi<double>()).margin(tolerance));
  CHECK( mx.variable_as_number("d") == Approx(boost::math::constants::pi<double>() / 180.0).margin(tolerance));
  CHECK( mx.variable_as_number("ee") == Approx(boost::math::constants::e<double>()).margin(tolerance));

  CHECK( mx.variable_as_number("f") == Approx(pconstants::me).margin(tolerance));
  CHECK( mx.variable_as_number("g") == Approx(pconstants::mp).margin(tolerance));
  CHECK( mx.variable_as_number("h") == Approx(pconstants::mmu).margin(tolerance));
  CHECK( mx.variable_as_number("i") == Approx(pconstants::c).margin(tolerance));
  CHECK( mx.variable_as_number("j") == Approx(pconstants::e).margin(tolerance));
}



TEST_CASE("newline_separation")
{
  string str = "x=1;\n y=2;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.variable_as_number("x") == 1 );
  CHECK( mx.variable_as_number("y") == 2 );
}

TEST_CASE("semicolon_separation")
{
  string str = "x=1;y=2;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.variable_as_number("x") == 1 );
  CHECK( mx.variable_as_number("y") == 2 );
}

TEST_CASE("floating_point")
{
  // 7/32 has an exact floating point representation
  string str = "x=1.234;y=0.21875; ze2=.21875e2; wep2=.21875E2; wep02=.21875e+02;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.variable_as_number("x") == Approx(1.234).margin(tolerance));
  CHECK( mx.variable_as_number("y") == 0.21875 );
  CHECK( mx.variable_as_number("ze2") == 21.875);
  CHECK( mx.variable_as_number("wep2") == 21.875);
  CHECK( mx.variable_as_number("wep02") == 21.875);
}

TEST_CASE("variable_assignment_expression")
{
  string str = "foo.bar=pi*sin(1.2e-4)^0.69;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.variable_as_number("foo.bar") ==
                     Approx(boost::math::constants::pi<double>() * pow(sin(1.2e-4), 0.69)).margin(tolerance));
}

TEST_CASE("caps_variable_assignment")
{
  string str = "X=1;Y=X;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.variable_as_number("y") == 1 );
}

// ==========================================================================
// commands

TEST_CASE("command")
{
  string str = "beam;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("beam");
  CHECK( cmd.name() == "beam" );
  CHECK( cmd.attribute_count() == 0 );
}

#if 0
// This is defined away because upper-case commands do not work.
TEST_CASE("upper_command")
{
  string str = "BEAM;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.command_count() == 1 );

  MadX_command cmd = mx.command(0);
  CHECK( cmd.name() == "beam" );
  CHECK( cmd.attribute_count() == 0 );
}
#endif

TEST_CASE("command_attrs")
{
  string str = "beam,a=1,B=3*(4+5);";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("beam");
  CHECK( cmd.name() == "beam" );
  CHECK( cmd.attribute_count() == 2 );
  CHECK( cmd.attribute_as_number("a") == 1 );
  CHECK( cmd.attribute_as_number("b") ==  3*(4+5) );
}

#if 0
TEST_CASE("command_str_attrs1")
{
  string str = "title, S = \"Tevatron Collider Run II Lattice\";";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.command_count() == 1 );

  MadX_command cmd = mx.command(0);
  CHECK( cmd.name() == "title" );
  CHECK( cmd.attribute_count() == 1 );
  CHECK( cmd.attribute_as_string("s") == "Tevatron Collider Run II Lattice" );
}
#endif

#if 0
TEST_CASE("command_str_attrs2")
{
  string str = "title, S = 'Tevatron Collider Run II Lattice';";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.command_count() == 1 );

  MadX_command cmd = mx.command(0);
  CHECK( cmd.name() == "title" );
  CHECK( cmd.attribute_count() == 1 );
  CHECK( cmd.attribute_as_string("s") == "Tevatron Collider Run II Lattice" );
}
#endif


TEST_CASE("command_particle_attrs")
{
  string str = "beAM, particle=proton;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("beam");
  CHECK( cmd.name() == "beam" );
  CHECK( cmd.attribute_count() == 5 );
  CHECK( cmd.attribute_as_string("particle") == "proton");
}

TEST_CASE("command_special_attrs1")
{
  string str = "mp: multipole, knl:={0, 1, 1}, type=octpn;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("mp");
  CHECK( cmd.name() == "multipole" );
  CHECK( cmd.attribute_count() == 2 );
  CHECK( cmd.attribute_as_string("type") == "octpn");
}

TEST_CASE("command_special_attrs2")
{
  string str = "mp: multipole, knl:={0, 1, 1}, TYPE=wgl;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("mp");
  CHECK( cmd.name() == "multipole" );
  CHECK( cmd.attribute_count() == 2 );
  CHECK( cmd.attribute_as_string("type") == "wgl");
}

TEST_CASE("command_special_attrs3")
{
  string str = "mp: multipole, knl:={0, 1, 1}, type=\"special\";";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("mp");
  CHECK( cmd.name() == "multipole" );
  CHECK( cmd.attribute_count() == 2 );
  CHECK( cmd.attribute_as_string("type") == "special");
}

TEST_CASE("command_omitted_comma")
{
  string str = "call file = './foo.dbx';";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.command_count() == 0 );

  CHECK( mx.variable_as_number("a") == 3 );
  CHECK( mx.variable_as_number("b") == 2 );
}

TEST_CASE("command_beam_particle")
{
  string str = " BEAM, PARTICLE=Proton, MASS=0.93827, CHARGE=1., ENERGY=0.93827 + 0.160;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );

  MadX_command cmd = mx.command("beam");
  CHECK( cmd.name() == "beam" );
  CHECK( cmd.attribute_count() == 6 );
  CHECK( cmd.attribute_as_number("mass") == Approx(0.93827));
  CHECK( cmd.attribute_as_number("charge") == 1 );
  CHECK( cmd.attribute_as_number("energy") == Approx(0.93827 + 0.160));
}

TEST_CASE("command_beam_particle_abbreviate")
{
  string str = " BEAM, PARTICLE=Prot, MASS=0.93827, CHARGE=1., ENERGY=0.93827 + 0.160;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );

  MadX_command cmd = mx.command("beam");
  CHECK( cmd.name() == "beam" );
  CHECK( cmd.attribute_count() == 6 );
  CHECK( cmd.attribute_as_number("mass") == Approx(0.93827));
  CHECK( cmd.attribute_as_number("charge") == 1 );
  CHECK( cmd.attribute_as_number("energy") == Approx(0.93827 + 0.160));
}

TEST_CASE("command_assign")
{
  string str = "q1: quadrupole,l=3.14;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("q1");
  CHECK( cmd.name() == "quadrupole" );
  CHECK( cmd.attribute_count() == 1 );
  CHECK( cmd.attribute_as_number("l") == Approx(3.14).margin(tolerance));
}

TEST_CASE("subscripted_ident")
{
  string str = "foo: quadrupole, a=1; x=foo->a;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.variable_as_number("x") == 1 );
}

TEST_CASE("vector_attribute")
{
    string str = "mpfi1: multipole,knl:={ 0, 1.1, 2.2, 3.3, 4.4, 5.5 };";
    MadX mx;
    CHECK_NOTHROW( parse_madx( str, mx ) );
    MadX_command cmd = mx.command("mpfi1");
    std::vector<double > knl(cmd.attribute_as_number_seq("knl"));

    CHECK(knl.size() == 6);
    const double tolerance = 1.0e-12;
    CHECK(knl.at(0) == Approx(0.0).margin(tolerance));
    CHECK(knl.at(1) == Approx(1.1).margin(tolerance));
    CHECK(knl.at(2) == Approx(2.2).margin(tolerance));
    CHECK(knl.at(3) == Approx(3.3).margin(tolerance));
    CHECK(knl.at(4) == Approx(4.4).margin(tolerance));
    CHECK(knl.at(5) == Approx(5.5).margin(tolerance));
}

TEST_CASE("continuation")
{
  string str = "q1: quadrupole,l=\n3.14,k1=0.2;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("q1");
  CHECK( cmd.name() == "quadrupole" );
  CHECK( cmd.attribute_count() ==  2 );
  CHECK( cmd.attribute_as_number("l") == Approx(3.14).margin(tolerance));
  CHECK( cmd.attribute_as_number("k1") == Approx(0.2).margin(tolerance));
}

TEST_CASE("continuation2")
{
  string str = "q2: quadrupole,l=! )junk)\n3.14,k1=0.2;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("q2");
  CHECK( cmd.name() ==  "quadrupole" );
  CHECK( cmd.attribute_count() == 2 );
  CHECK( cmd.attribute_as_number("l") == Approx(3.14).margin(tolerance));
  CHECK( cmd.attribute_as_number("k1") == Approx(0.2).margin(tolerance));
}

TEST_CASE("continuation3")
{
  string str = "q3: quadrupole,l=3.14,\nk1=0.2;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("q3");
  CHECK( cmd.name() == "quadrupole" );
  CHECK( cmd.attribute_count() == 2 );
  CHECK( cmd.attribute_as_number("l") == Approx(3.14).margin(tolerance));
  CHECK( cmd.attribute_as_number("k1") == Approx(0.2).margin(tolerance));
}

TEST_CASE("tkicker")
{
  string str = "tk: tkicker, hkick=0.01;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("tk");
  CHECK( cmd.name() == "tkicker" );
  CHECK( cmd.attribute_count() ==  1 );
  CHECK( cmd.attribute_as_number("hkick") == Approx(0.01).margin(tolerance));
}

TEST_CASE("constfoc")
{
  string str = "cf: constfoc, betaH=0.01, betaV=0.02, betaL=0.03, nuL=0.04;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("cf");
  CHECK( cmd.name() == "constfoc" );
  CHECK( cmd.attribute_count() == 4 );
  CHECK( cmd.attribute_as_number("betaH") == Approx(0.01).margin(tolerance));
  CHECK( cmd.attribute_as_number("betaV") == Approx(0.02).margin(tolerance));
  CHECK( cmd.attribute_as_number("betaL") == Approx(0.03).margin(tolerance));
  CHECK( cmd.attribute_as_number("nuL") == Approx(0.04).margin(tolerance));
}

TEST_CASE("comment_at_tail")
{
  string str = "tk: tkicker, hkick=0.01;!afdsak";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
}

TEST_CASE("parse_file")
{
  MadX   mx;

  //parse_madx_file( "./PS.madx", mx );
  //CHECK_NOTHROW( parse_madx_file( "./pstest/PS.madx", mx ) );
}


TEST_CASE("matrix")
{
  string str = "m1: matrix, type=abc, L=1.2, kick1=0.001, kick2=0.001, kick3=0.002, kick4=0.002, kick5=0.003, kick6=0.003, rm11=1.1, rm12=1.2, rm32=3.2, rm54=5.4, tm111=1.11, tm321=3.21;";

  MadX mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.label_count() == 1 );

  MadX_command cmd = mx.command("m1");
  CHECK( cmd.name() == "matrix" );
  CHECK( cmd.attribute_count() == 14 );
  CHECK( cmd.attribute_as_string("type") == "abc" );
  CHECK( cmd.attribute_as_number("l") == Approx(1.2).margin(tolerance));
  CHECK( cmd.attribute_as_number("kick1") == Approx(0.001).margin(tolerance));
  CHECK( cmd.attribute_as_number("kick2") == Approx(0.001).margin(tolerance));
  CHECK( cmd.attribute_as_number("kick3") == Approx(0.002).margin(tolerance));
  CHECK( cmd.attribute_as_number("kick4") == Approx(0.002).margin(tolerance));
  CHECK( cmd.attribute_as_number("kick5") == Approx(0.003).margin(tolerance));
  CHECK( cmd.attribute_as_number("kick6") == Approx(0.003).margin(tolerance));
  CHECK( cmd.attribute_as_number("rm11") == Approx(1.1).margin(tolerance));
  CHECK( cmd.attribute_as_number("rm12") == Approx(1.2).margin(tolerance));
  CHECK( cmd.attribute_as_number("rm32") == Approx(3.2).margin(tolerance));
  CHECK( cmd.attribute_as_number("rm54") == Approx(5.4).margin(tolerance));
  CHECK( cmd.attribute_as_number("tm111") == Approx(1.11).margin(tolerance));
  CHECK( cmd.attribute_as_number("tm321") == Approx( 3.21).margin(tolerance));
}

#if 0
TEST_CASE("continuation4")
{
  string str = "twiss, save,   betx=28.871,alfx=-0.069,mux=0.0,dx=2.682,dpx=-0.073,\nbety= 5.264,alfy=-0.006,muy=0.0,dy=0.0,dpy=0.0;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
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
#endif

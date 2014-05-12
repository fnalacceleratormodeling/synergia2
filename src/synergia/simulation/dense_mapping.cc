#include "dense_mapping.h"

Dense_mapping::Dense_mapping() :
    constant(boost::extents[6]), linear(boost::extents[6][6])
{
}

Dense_mapping::Dense_mapping(Fast_mapping const& fast_mapping) :
    constant(boost::extents[6]), linear(boost::extents[6][6])
{
    for (int i = 0; i < 6; ++i) {
        constant[i] = 0.0;
        for (int j = 0; j < 6; ++j) {
            linear[i][j] = 0.0;
        }
    }
    std::vector<std::vector<Fast_mapping_terms > > const& terms =
            fast_mapping.get_terms();
    for (int i = 0; i < 6; ++i) {
        Fast_mapping_terms::const_iterator telem;
        for (telem = terms[i][0].begin(); telem != terms[i][0].end(); ++telem) {
            constant[i] += telem->coeff();
        }
        for (telem = terms[i][1].begin(); telem != terms[i][1].end(); ++telem) {
            linear[i][telem->index(0)] = telem->coeff();
        }
    }
    length = fast_mapping.get_length();
}

double
Dense_mapping::get_length() const
{
    return length;
}

MArray1d_ref
Dense_mapping::get_constant_term() const
{
    return constant;
}

MArray2d_ref
Dense_mapping::get_linear_term() const
{
    return linear;
}

MArray2d
Dense_mapping::get_linear_term_mad()
{
    const int magic_number = 3;
    Jet__environment::BeginEnvironment(magic_number);
    coord x(0.0);
    coord xp(0.0);
    coord y(0.0);
    coord yp(0.0);
    coord cdt(0.0);
    coord dpop(0.0);
    Jet__environment::EndEnvironment();
  
#if 0  
    // Phi is the transformation from the
    // "old" to the "new" coordinates.
    // ----------------------------------
    Mapping Phi;
    Phi.SetComponent( 0, x + x*x );  
    Phi.SetComponent( 1, y + y*y );  
  
    // T is a transformation expressed in
    // terms of the "old" coordinates.
    // Here, it is merely a scale transformation.
    // ------------------------------------------
    Mapping T;
    T.SetComponent( 0, 2.0*x );
    T.SetComponent( 1, 2.0*y );
  
  
    // We create a second environment at the
    // image of the reference point under the
    // transformation T.
    // --------------------------------------
    Vector new_reference( T( Phi.getReference() ) );
  
    Jet__environment::BeginEnvironment(9);
    coord xx(new_reference[0]);
    coord yy(new_reference[1]);
    Jet__environment::EndEnvironment();
  
  
    // Phi_2 is the same as Phi - the coordinate
    // transformation from "old" to "new"
    // coordinates - but viewed from the T-imaged
    // reference point.
    // 
    // Because the environment was just created,
    // it is currently the default environment,
    // so need not be explicitly passed to
    // constructors.
    // ------------------------------------------
    Mapping Phi_2;
    Phi_2.SetComponent( 0, xx + xx*xx );
    Phi_2.SetComponent( 1, yy + yy*yy );
  
  
    // Here, I'll store the inverse of
    // Phi_2 in a variable, though this is
    // not necessary. (see on_origin code)
    // -----------------------------------
    Mapping Psi( Phi_2.inverse() );
  
  
    // F is the transformation T
    // expressed in the "new" coordinates
    // ----------------------------------
    Mapping F(   Phi( T( Psi ) )   );
#endif  
  
    return linear;
}

Dense_mapping::~Dense_mapping()
{

}


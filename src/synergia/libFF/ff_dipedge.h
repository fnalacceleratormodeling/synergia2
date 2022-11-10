#ifndef FF_DIPEDGE_H
#define FF_DIPEDGE_H

include "synergia/libFF/ff_algorithm.h"
#include "synergia/utils/simple_timer.h"

class FF_dipedge : public FF_element
{
private:
    double get_reference_cdt(double re_2_1, double re_4_3, double * te, Reference_particle & reference_particle);
public:
    FF_dipedge();
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_dipedge();
};

typedef boost::shared_ptr<FF_dipedge > FF_dipedge_sptr;


#
// =======================================================================================
// ======================  END OF OLD DIPEDGE CODE =======================================
// =======================================================================================


// p [Gev/c] = -- * B*rho [ Tesla meters ]
#define PH_CNV_brho_to_p   (1.0e-9 * pconstants::c)

// Calculations of dipedge taken from MAD-X SUBROUTINE tmfrng from
// file twiss.f90.

// Comments from that file:
/*!----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for fringe field of a dipole.                      *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     h         (double)  curvature of magnet body.                    *
  !     sk1       (double)  quadrupole strength in magnet body.          *
  !     edge      (double)  edge focussing angle.                        *
  !     he        (double)  curvature of pole face.                      *
  !     sig       (double)  sign: +1 for entry, -1 for exit.             *
  !     corr      (double)  correction factor according to slac 75.      *
  !     Output:                                                          *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second order terms.                          *
*/

/* MAD-X Documentation on dipedge:

label:
DIPEDGE, H=real, E1=real, FINT=real,
HGAP=real, TILT=real;

A DIPEDGE  has zero length and five attributes.
H  Is angle/length or 1/rho (default: 0 m−1 - for the default the dipedge element
has no effect). (must be equal to that of the associated SBEND)

E1 The rotation angle for the pole face. The sign convention is as for a SBEND
bending magnet. Note that it is different for an entrance and an exit. (default:
0 rad).

FINT  field integral as for SBEND sector bend. Note that each DIPEDGE has its own
FINT, so that specifying FINTX is no longer necessary.

HGAP  half gap height of the associated SBEND bending magnet.

TILT  The roll angle about the longitudinal axis (default: 0 rad, i.e. a horizontal
bend). A positive angle represents a clockwise rotation.
 */



namespace dipedge_impl
{
    struct DipedgeParams
    {

        double e1;
        double h;
        double fint;
        double hgap;
        double tilt;

        double scale;

        double re_2_1;
        double re_4_3;
        double te[10];

        double strength;

        double pref_b;
        double m_b;
        double ref_cdt;

        double ce1, se1;
        double tanedge;
        double secedge;
        double sn;
        double psip;
        double re_2_1;
        double re_4_3;

    };

    template<class BP> // BP is bunch particles?
    struct PropDipedge
    {
        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        const DipedgeParams dipedge_params;

        PropDipedge(BP & bp, DipedgeParams const& dipedge_params)
            : p(bp.parts), masks(bp.masks), dipedge_params(dipedge_params)
            { }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (masks(i))
            {
                if (sp.ledge)
                {
                    // slot
                    FF_algorithm::slot_unit( 
                            p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                            sp.ce1, sp.se1, sp.pref_b, sp.m_b );

                    // edge:
                    // 1. chef fixed angle (only kicks yp)
                    // FF_algorithm::edge_unit(
                    //        p(i,2), p(i,3), us_edge_k );
                    //
                    // 2. chef per-particle angle
                    // FF_algorithm::edge_unit(
                    //        p(i,2), p(i,1), p(i,3), dpop, us_edge_k_p);
                    //
                    // 3. ref particle angle (kicks both xp and yp)
                    FF_algorithm::edge_unit(
                            p(i,2), p(i,1), p(i,3), 
                            sp.us_edge_k_x, sp.us_edge_k_y, 0); 
                }

                // bend
                FF_algorithm::bend_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                        dphi, sp.strength, sp.pref_b, sp.m_b, sp.ref_cdt, 
                        phase, term );

                if (sp.redge)
                {
                    // edge
                    // FF_algorithm::edge_unit(y, yp, ds_edge_k);
                    // FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);
                    FF_algorithm::edge_unit(
                            p(i,2), p(i,1), p(i,3), 
                            sp.ds_edge_k_x, sp.ds_edge_k_y, 0); 

                    // slot
                    FF_algorithm::slot_unit(
                            p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                            sp.ce2, sp.se2, sp.pref_b, sp.m_b);
                }
            }
        }
    };

    template<class BP>
    struct PropDipedgeSimd
    {
        using gsv_t = typename BP::gsv_t;
        using parts_t = typename BP::parts_t;

        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        const SbendParams sp;

        const double dphi;
        const Kokkos::complex<double> phase;
        const Kokkos::complex<double> term;

        PropSbendSimd(BP & bp, SbendParams const& sp)
            : p(bp.parts), masks(bp.masks), sp(sp)
            , dphi( -(sp.angle - (sp.e1 + sp.e2)) )   // -psi
            , phase(std::exp(std::complex<double>(0.0, -dphi)))
            , term( std::complex<double>(0.0, sp.length / sp.angle) *
                    std::complex<double>(1.0 - cos(sp.angle), - sin(sp.angle)) *
                    std::complex<double>(cos(sp.e2), -sin(sp.e2)) )
        { }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int idx) const
        {
            int i = idx * gsv_t::size();

            int m = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) m |= masks(x);

            if (m)
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                if (sp.ledge)
                {
                    // slot
                    FF_algorithm::slot_unit<gsv_t>( 
                            p0, p1, p2, p3, p4, p5,
                            sp.ce1, sp.se1, sp.pref_b, sp.m_b );

                    // edge:
                    // 1. chef fixed angle (only kicks yp)
                    // FF_algorithm::edge_unit<gsv_t>(
                    //        p2, p3, sp.us_edge_k );
                    //
                    // 2. chef per-particle angle
                    // FF_algorithm::edge_unit<gsv_t>(
                    //        p2, p1, p3, p5, sp.us_edge_k_p);
                    //
                    // 3. ref particle angle (kicks both xp and yp)
                    FF_algorithm::edge_unit<gsv_t>( 
                            p2, p1, p3, sp.us_edge_k_x, sp.us_edge_k_y, 0); 
                }

                // bend
                FF_algorithm::bend_unit<gsv_t>(
                        p0, p1, p2, p3, p4, p5,
                        dphi, sp.strength, sp.pref_b, sp.m_b, sp.ref_cdt, 
                        phase, term );

                if (sp.redge)
                {
                    // edge
                    // FF_algorithm::edge_unit<gsv_t>(p2, p3, sp.ds_edge_k);
                    // FF_algorithm::edge_unit<gsv_t>(p2, p1, p3, p5, sp.ds_edge_k_p);
                    FF_algorithm::edge_unit<gsv_t>( 
                            p2, p1, p3, sp.ds_edge_k_x, sp.ds_edge_k_y, 0); 

                    // slot
                    FF_algorithm::slot_unit<gsv_t>(
                            p0, p1, p2, p3, p4, p5,
                            sp.ce2, sp.se2, sp.pref_b, sp.m_b);
                }

                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
                p5.store(&p(i, 5));
            }
        }
    };

    namespace fa  = ::FF_algorithm;

    template<class BP>
    struct PropSbendCF
    {
        using gsv_t = typename BP::gsv_t;
        using parts_t = typename BP::parts_t;
        using const_masks_t = typename BP::const_masks_t;

        parts_t p;
        const_masks_t masks;
        const SbendParams sp;

        const double step_angle;
        const double step_ref_cdt;

        const double step_kl[10];

        const Kokkos::complex<double> phase_e1;
        const Kokkos::complex<double> phase_e2;

        const Kokkos::complex<double> step_phase[4];
        const Kokkos::complex<double> step_term[4];
        const double step_dphi[4];

        PropSbendCF( BP & bp, SbendParams const& sp )
            : p(bp.parts), masks(bp.masks), sp(sp)

            , step_angle(sp.angle/sp.steps)
            , step_ref_cdt(sp.ref_cdt/sp.steps)

            , step_kl{  // k_b[i] = k_l[i] * scale
                sp.k_l[0] * sp.scale * sp.length / sp.steps,
                sp.k_l[1] * sp.scale * sp.length / sp.steps,
                sp.k_l[2] * sp.scale * sp.length / sp.steps,
                sp.k_l[3] * sp.scale * sp.length / sp.steps,
                sp.k_l[4] * sp.scale * sp.length / sp.steps,
                sp.k_l[5] * sp.scale * sp.length / sp.steps,
                sp.k_l[6] * sp.scale * sp.length / sp.steps,
                sp.k_l[7] * sp.scale * sp.length / sp.steps,
                sp.k_l[8] * sp.scale * sp.length / sp.steps,
                sp.k_l[9] * sp.scale * sp.length / sp.steps }

            , phase_e1(fa::bend_edge_phase(sp.e1))
            , phase_e2(fa::bend_edge_phase(sp.e2))

            , step_phase{
                fa::sbend_unit_phase(by6::c1, step_angle),
                fa::sbend_unit_phase(by6::c2, step_angle),
                fa::sbend_unit_phase(by6::c3, step_angle),
                fa::sbend_unit_phase(by6::c4, step_angle) }

            , step_term{
                fa::sbend_unit_term(by6::c1, step_angle, sp.r0),
                fa::sbend_unit_term(by6::c2, step_angle, sp.r0),
                fa::sbend_unit_term(by6::c3, step_angle, sp.r0),
                fa::sbend_unit_term(by6::c4, step_angle, sp.r0) }

            , step_dphi{
                fa::sbend_dphi(by6::c1, step_angle),
                fa::sbend_dphi(by6::c2, step_angle),
                fa::sbend_dphi(by6::c3, step_angle),
                fa::sbend_dphi(by6::c4, step_angle) }
        { }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (masks(i))
            {
                if (sp.ledge)
                {
                    // slot
                    FF_algorithm::slot_unit(
                            p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                            sp.ce1, sp.se1, sp.pref_b, sp.m_b );

                    // edge
                    //FF_algorithm::edge_unit(y, yp, us_edge_k);
                    //FF_algorithm::edge_unit(y, xp, yp, dpop, us_edge_k_p);
                    FF_algorithm::edge_unit(
                            p(i,2), p(i,1), p(i,3), 
                            sp.us_edge_k_x, sp.us_edge_k_y, 0); 

                    // bend edge (thin, but with face angle)
                    FF_algorithm::bend_edge(
                            p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                            sp.e1, phase_e1, sp.strength, sp.pref_b, sp.m_b);
                }

                // bend body
                FF_algorithm::bend_yoshida6<typename BP::part_t,
                    fa::thin_cf_kick_5<typename BP::part_t>, 5> ( 
                            p(i,0), p(i,1), p(i,2), 
                            p(i,3), p(i,4), p(i,5),
                            sp.pref_b, sp.m_b, step_ref_cdt,
                            step_kl, step_dphi, step_phase, step_term,
                            sp.r0, sp.strength, sp.steps );

                if (sp.redge)
                {
                    // bend edge (thin, but with face angle)
                    FF_algorithm::bend_edge(
                            p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                            sp.e2, phase_e2, sp.strength, sp.pref_b, sp.m_b);

                    // edge
                    //FF_algorithm::edge_unit(y, yp, ds_edge_k);
                    //FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);
                    FF_algorithm::edge_unit(
                            p(i,2), p(i,1), p(i,3), 
                            sp.ds_edge_k_x, sp.ds_edge_k_y, 0); 

                    // slot
                    FF_algorithm::slot_unit(
                            p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                            sp.ce2, sp.se2, sp.pref_b, sp.m_b);
                }
            }
        }
    };

    template<class BP>
    struct PropSbendCFSimd
    {
        using gsv_t = typename BP::gsv_t;
        using parts_t = typename BP::parts_t;
        using const_masks_t = typename BP::const_masks_t;

        parts_t p;
        const_masks_t masks;
        const SbendParams sp;

        const double step_angle;
        const double step_ref_cdt;

        const double step_kl[10];

        const Kokkos::complex<double> phase_e1;
        const Kokkos::complex<double> phase_e2;

        const Kokkos::complex<double> step_phase[4];
        const Kokkos::complex<double> step_term[4];
        const double step_dphi[4];

        PropSbendCFSimd( BP & bp, SbendParams const& sp )
            : p(bp.parts), masks(bp.masks), sp(sp)

            , step_angle(sp.angle/sp.steps)
            , step_ref_cdt(sp.ref_cdt/sp.steps)

            , step_kl{  // k_b[i] = k_l[i] * scale
                sp.k_l[0] * sp.scale * sp.length / sp.steps,
                sp.k_l[1] * sp.scale * sp.length / sp.steps,
                sp.k_l[2] * sp.scale * sp.length / sp.steps,
                sp.k_l[3] * sp.scale * sp.length / sp.steps,
                sp.k_l[4] * sp.scale * sp.length / sp.steps,
                sp.k_l[5] * sp.scale * sp.length / sp.steps,
                sp.k_l[6] * sp.scale * sp.length / sp.steps,
                sp.k_l[7] * sp.scale * sp.length / sp.steps,
                sp.k_l[8] * sp.scale * sp.length / sp.steps,
                sp.k_l[9] * sp.scale * sp.length / sp.steps }

            , phase_e1(fa::bend_edge_phase(sp.e1))
            , phase_e2(fa::bend_edge_phase(sp.e2))

            , step_phase{
                fa::sbend_unit_phase(by6::c1, step_angle),
                fa::sbend_unit_phase(by6::c2, step_angle),
                fa::sbend_unit_phase(by6::c3, step_angle),
                fa::sbend_unit_phase(by6::c4, step_angle) }

            , step_term{
                fa::sbend_unit_term(by6::c1, step_angle, sp.r0),
                fa::sbend_unit_term(by6::c2, step_angle, sp.r0),
                fa::sbend_unit_term(by6::c3, step_angle, sp.r0),
                fa::sbend_unit_term(by6::c4, step_angle, sp.r0) }

            , step_dphi{
                fa::sbend_dphi(by6::c1, step_angle),
                fa::sbend_dphi(by6::c2, step_angle),
                fa::sbend_dphi(by6::c3, step_angle),
                fa::sbend_dphi(by6::c4, step_angle) }
        { }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int idx) const
        {
            int i = idx * gsv_t::size();

            int m = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) m |= masks(x);

            if (m)
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                if (sp.ledge)
                {
                    // slot
                    FF_algorithm::slot_unit<gsv_t>(
                            p0, p1, p2, p3, p4, p5,
                            sp.ce1, sp.se1, sp.pref_b, sp.m_b );

                    // edge
                    //FF_algorithm::edge_unit(y, yp, us_edge_k);
                    //FF_algorithm::edge_unit(y, xp, yp, dpop, us_edge_k_p);
                    FF_algorithm::edge_unit<gsv_t>(
                            p2, p1, p3, 
                            sp.us_edge_k_x, sp.us_edge_k_y, 0); 

                    // bend edge (thin, but with face angle)
                    FF_algorithm::bend_edge<gsv_t>(
                            p0, p1, p2, p3, p4, p5,
                            sp.e1, phase_e1, sp.strength, sp.pref_b, sp.m_b);
                }

                // bend body
                FF_algorithm::bend_yoshida6<gsv_t, 
                    fa::thin_cf_kick_5<gsv_t>, 5> ( 
                            p0, p1, p2, p3, p4, p5,
                            sp.pref_b, sp.m_b, step_ref_cdt,
                            step_kl, step_dphi, step_phase, step_term,
                            sp.r0, sp.strength, sp.steps );

                if (sp.redge)
                {
                    // bend edge (thin, but with face angle)
                    FF_algorithm::bend_edge<gsv_t>(
                            p0, p1, p2, p3, p4, p5,
                            sp.e2, phase_e2, sp.strength, sp.pref_b, sp.m_b);

                    // edge
                    //FF_algorithm::edge_unit(y, yp, ds_edge_k);
                    //FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);
                    FF_algorithm::edge_unit<gsv_t>(
                            p2, p1, p3, 
                            sp.ds_edge_k_x, sp.ds_edge_k_y, 0); 

                    // slot
                    FF_algorithm::slot_unit<gsv_t>(
                            p0, p1, p2, p3, p4, p5,
                            sp.ce2, sp.se2, sp.pref_b, sp.m_b);
                }

                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
                p5.store(&p(i, 5));
            }
        }
    };

    inline void prop_reference(
            Reference_particle & ref_l, 
            DipedgeParams & dipedge_params )
    {
        double pref_l = ref_l.get_momentum();
        double    m_l = ref_l.get_mass();

        // propagate the reference particle, and set the edge kick strength 
        // from the reference particle
        double    x_l = ref_l.get_state()[Bunch::x];
        double   xp_l = ref_l.get_state()[Bunch::xp];
        double    y_l = ref_l.get_state()[Bunch::y];
        double   yp_l = ref_l.get_state()[Bunch::yp];
        double  cdt_l = 0.0;
        double dpop_l = ref_l.get_state()[Bunch::dpop];

        double dphi =  -(sp.angle - (sp.e1 + sp.e2));

        Kokkos::complex<double> phase = 
            std::exp(std::complex<double>(0.0, -dphi));

        Kokkos::complex<double> term  = 
            std::complex<double>(0.0, sp.length / sp.angle) *
            std::complex<double>(1.0 - cos(sp.angle), - sin(sp.angle)) *
            std::complex<double>(cos(sp.e2), -sin(sp.e2));

        if (sp.ledge)
        {
            // slot
            FF_algorithm::slot_unit(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.ce1, sp.se1, pref_l, m_l);

            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            sp.us_edge_k_x = sp.us_edge_k_p * (xp_l/zp_l);
            sp.us_edge_k_y = sp.us_edge_k_p * (yp_l/zp_l);

            // edge kick that is in accordance with 
            // method 1 (chef fixed angle)
            // double brho_l = pref_l / PH_CNV_brho_to_p;
            // sp.us_edge_k = sp.strength * tan(atan2(xp_l, zp_l)) / brho_l;
            // FF_algorithm::edge_unit(y_l, yp_l, sp.us_edge_k );

            // edge kick strenth are scaled to bunch. so need to div by "scale" to scale
            // it to the lattice reference
            // in accordance with method 3
            FF_algorithm::edge_unit(y_l, xp_l, yp_l, 
                    sp.us_edge_k_x/sp.scale, sp.us_edge_k_y/sp.scale, 0);

        }

        FF_algorithm::bend_complete(
                x_l, xp_l, y_l, yp_l, cdt_l, dpop_l,
                dphi, sp.strength, pref_l, m_l, 0.0/*ref cdt*/, 
                phase, term);

        if (sp.redge)
        {
            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            sp.ds_edge_k_x = sp.ds_edge_k_p * (xp_l/zp_l);
            sp.ds_edge_k_y = sp.ds_edge_k_p * (yp_l/zp_l);

            // edge kick that is in accordance with 
            // method 1 (chef fixed angle)
            // double brho_l = pref_l / PH_CNV_brho_to_p;
            // sp.ds_edge_k = - sp.strength * tan(atan2(xp_l, zp_l)) / brho_l;
            // FF_algorithm::edge_unit(y_l, yp_l, sp.ds_edge_k );

            // edge kick strenth are scaled to bunch. so need to div by "scale" to scale
            // it to the lattice reference
            FF_algorithm::edge_unit(y_l, xp_l, yp_l, 
                    sp.ds_edge_k_x/sp.scale, sp.ds_edge_k_y/sp.scale, 0);

            // slot
            FF_algorithm::slot_unit(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.ce2, sp.se2, pref_l, m_l);
        }

        ref_l.set_state(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l);
        sp.ref_cdt = cdt_l;
    }

    inline void prop_reference_cf(
            Reference_particle & ref_l, 
            SbendParams & sp )
    {
        double pref_l = ref_l.get_momentum();
        double    m_l = ref_l.get_mass();

        // propagate the reference particle, and set the edge kick strength 
        // from the reference particle
        double    x_l = ref_l.get_state()[Bunch::x];
        double   xp_l = ref_l.get_state()[Bunch::xp];
        double    y_l = ref_l.get_state()[Bunch::y];
        double   yp_l = ref_l.get_state()[Bunch::yp];
        double  cdt_l = 0.0;
        double dpop_l = ref_l.get_state()[Bunch::dpop];

        double step_length = sp.length / sp.steps;
        double step_angle = sp.angle/sp.steps;

        // step strength for reference
        double step_kl[10];
        for (int i=0; i<10; ++i) step_kl[i] = sp.k_l[i] * step_length;

        Kokkos::complex<double> phase_e1 = FF_algorithm::bend_edge_phase(sp.e1);
        Kokkos::complex<double> phase_e2 = FF_algorithm::bend_edge_phase(sp.e2);

        Kokkos::complex<double> step_phase[4] = { 
            fa::sbend_unit_phase(by6::c1, step_angle),
            fa::sbend_unit_phase(by6::c2, step_angle),
            fa::sbend_unit_phase(by6::c3, step_angle),
            fa::sbend_unit_phase(by6::c4, step_angle) 
        };

        Kokkos::complex<double> step_term[4] = {
            fa::sbend_unit_term(by6::c1, step_angle, sp.r0),
            fa::sbend_unit_term(by6::c2, step_angle, sp.r0),
            fa::sbend_unit_term(by6::c3, step_angle, sp.r0),
            fa::sbend_unit_term(by6::c4, step_angle, sp.r0)
        };

        double step_dphi[4] = {
            fa::sbend_dphi(by6::c1, step_angle),
            fa::sbend_dphi(by6::c2, step_angle),
            fa::sbend_dphi(by6::c3, step_angle),
            fa::sbend_dphi(by6::c4, step_angle),
        };

        if (sp.ledge)
        {
            // slot
            FF_algorithm::slot_unit(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.ce1, sp.se1, pref_l, m_l);

            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            sp.us_edge_k_x = sp.us_edge_k_p * (xp_l/zp_l);
            sp.us_edge_k_y = sp.us_edge_k_p * (yp_l/zp_l);

            // edge
            //FF_algorithm::edge_unit(y_l, yp_l, us_edge_k/scale);
            //FF_algorithm::edge_unit(y_l, xp_l, yp_l, dpop_l, us_edge_k_p/scale);
            FF_algorithm::edge_unit(y_l, xp_l, yp_l, 
                    sp.us_edge_k_x/sp.scale, sp.us_edge_k_y/sp.scale, 0);

            // bend edge (thin)
            FF_algorithm::bend_edge(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.e1, phase_e1, sp.strength, pref_l, m_l);
        }

        FF_algorithm::bend_yoshida6<double, 
            FF_algorithm::thin_cf_kick_5<double>, 5> ( 
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l,
                    pref_l, m_l, 0.0 /* step ref_cdt */,
                    step_kl, step_dphi, step_phase, step_term,
                    sp.r0, sp.strength, sp.steps);

        if (sp.redge)
        {
            // bend edge (thin)
            FF_algorithm::bend_edge(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.e2, phase_e2, sp.strength, pref_l, m_l);

            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            sp.ds_edge_k_x = sp.ds_edge_k_p * (xp_l/zp_l);
            sp.ds_edge_k_y = sp.ds_edge_k_p * (yp_l/zp_l);

            // edge
            //FF_algorithm::edge_unit(y_l, yp_l, ds_edge_k);
            //FF_algorithm::edge_unit(y_l, xp_l, yp_l, dpop_l, ds_edge_k_p/scale);
            FF_algorithm::edge_unit(y_l, xp_l, yp_l, 
                    sp.ds_edge_k_x/sp.scale, sp.ds_edge_k_y/sp.scale, 0);

            // slot
            FF_algorithm::slot_unit(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.ce2, sp.se2, pref_l, m_l);
        }

        ref_l.set_state(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l);
        sp.ref_cdt = cdt_l;
    }

    // max multipole moments
    constexpr const int max_mp_order = 5;

    template<class BP>
    struct PropThinPoleSimd
    {
        using gsv_t = typename BP::gsv_t;

        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        kt::arr_t<double, 2*max_mp_order> k;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int idx) const
        {
            int i = idx * gsv_t::size();

            int m = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) m |= masks(x);

            if (m)
            {
                gsv_t  x(&p(i, 0));
                gsv_t xp(&p(i, 1));
                gsv_t  y(&p(i, 2));
                gsv_t yp(&p(i, 3));

                if (k[2] || k[3])   
                    FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, k.data+2);

                if (k[4] || k[5])   
                    FF_algorithm::thin_sextupole_unit(x, xp, y, yp, k.data+4);

                if (k[6] || k[7])   
                    FF_algorithm::thin_octupole_unit(x, xp, y, yp, k.data+6);

                if (k[8] || k[9])   
                    FF_algorithm::thin_magnet_unit(x, xp, y, yp, k.data+8, 5);

#if 0
                if (k[10] || k[11]) 
                    FF_algorithm::thin_magnet_unit(x, xp, y, yp, k.data+10, 6);

                if (k[12] || k[13]) 
                    FF_algorithm::thin_magnet_unit(x, xp, y, yp, k.data+12, 7);

                if (k[14] || k[15]) 
                    FF_algorithm::thin_magnet_unit(x, xp, y, yp, k.data+14, 8);
#endif

                xp.store(&p(i, 1));
                yp.store(&p(i, 3));
            }
        }
    };


}

namespace FF_dipedge
{

template<class BunchT>
inline void apply(Lattice_element_slice const& slice, BunchT & bunch)
{
    using namespace dipedge_impl;

    scoped_simple_timer timer("libFF_dipedge");

    auto const& ele = slice.get_lattice_element();

    DipedgeParams dipedge_params;

    dipedge_params.e1 = ele.get_double_attribute("e1", 0.0);
    dipedge_params.h = ele.get_double_attribute("h", 0.0);
    dipedge_params.fint = ele.get_double_attribute("fint", 0.0);
    dipedge_params.hgap = ele.get_double_attribute("hgap", 0.0);
    dipedge_params.tilt = ele.get_double_attribute("tilt", 0.0);
    dipedge_params.tanedg = tan(dipedge_params.e1);
    dipedge_params.secedg = 1.0 / cos(dipedge_params.e1);
    dipedge_params.sn = sin(dipedge_params.e1);
    dipedge_params.psip = 
    double tanedg =   tan(edge);
    double secedg =   1.0 / cos(edge);
    double sn     =   sin(edge);
    double psip   =   edge - corr * secedg * (1.0 + sn * sn);
    double re_2_1 =   h * tanedg;
    double re_4_3 = - h * tan(psip);
    Reference_particle       & ref_l = bunch.get_design_reference_particle();
    Reference_particle const & ref_b = bunch.get_reference_particle();

     dipedge_params.scale = ref_l.get_momentum() / 
        (ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]));


    // lattice reference
    double pref_l = ref_l.get_momentum();
    double brho_l = pref_l / PH_CNV_brho_to_p;
    double    m_l = ref_l.get_mass();
    int  charge_l = ref_l.get_charge();

    // the particle dpop is with respect to this momentum which goes with the bunch
    double pref_b = ref_b.get_momentum();
    double brho_b = pref_b / PH_CNV_brho_to_p;
    double    m_b = bunch.get_mass();

    // common
    dipedge_params.strength = brho_l * a / l;
    dipedge_params.pref_b   = pref_b;
    dipedge_params.m_b      = m_b;

    dipedge_params.ce1 = cos(-dipedge_params.e1);
    dipedge_params.se1 = sin(-dipedge_params.e1);

    using namespace Kokkos;
    using exec = typename BunchT::exec_space;

    auto apply = [&](ParticleGroup pg) {
        auto bp = bunch.get_bunch_particles(pg);
        if (!bp.num_valid()) return;

        // dipedge is a 0 length element with no multipole moments (currently)
#if LIBFF_USE_GSV
        prop_reference(ref_l, sp);

        auto range = RangePolicy<exec>(0, bp.size_in_gsv());
        PropDipedgeSimd<typename BunchT::bp_t> dipedge(bp, dipedge_params);
        Kokkos::parallel_for(range, dipedge);

#else // LIBFF_USE_GSV
        auto range = RangePolicy<exec>(0, bp.size());
        PropDipedge<typename BunchT::bp_t> dipedge(bp, dipedge_params);
        Kokkos::parallel_for(range, dipedge);
#endif // LIBFF_USE_GSV
    };

    apply(ParticleGroup::regular);
    apply(ParticleGroup::spectator);
    }

    // dipedge has 0 length so there is no need to increment trajectory here
    // bunch.get_reference_particle().increment_trajectory(0);

    Kokkos::fence();
}

}

#endif // FF_DIPEDGE_H


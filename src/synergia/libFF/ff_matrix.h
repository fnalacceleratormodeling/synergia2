#ifndef FF_MATRIX_H
#define FF_MATRIX_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/utils/simple_timer.h"


// Implementation of MAD-X arbitrary matrix element

/* MAD-X Documentation on matrix:

Arbitrary Matrix Element

label: MATRIX, TYPE=string, L=real,
    KICK1=real,..., KICK6=real,
    RM11=real, ..., RM66=real,
   TM111=real, ..., TM666=real;

The MATRIX element allows the definition of an arbitrary transfer matrix. It has four real
array attributes:
LLength of the element, which may be zero.

KICKi  Defines the kick of the element acting on the six phase space coordinates.

RMik   Defines the linear transfer matrix (6 × 6 terms) of the element.

TMikl  Defines the second-order terms (6 × 6 × 6 terms) of the element.

Data values that are not explicitly entered are taken from the identity transformation for the
RMik matrix elements, and taken as zero for the KICKi kick factors and the TMikl second order
terms. In the thin-lens tracking module a non-zero length for an arbitrary matrix is accepted,98
CHAPTER 11. ELEMENT TYPES
however no non-zero second order terms are allowed to avoid non symplectic tracking runs.
In the latter case the tracking run is aborted.label:

 */



// p [Gev/c] = -- * B*rho [ Tesla meters ]
#define PH_CNV_brho_to_p   (1.0e-9 * pconstants::c)

namespace matrix_impl
{
    struct MatrixParams
    {

        double  k1, k2, k3, k4, k5, k6; // KICKn parameters
        double  rm11, rm12, rm13, rm14, rm15, rm16,
                rm21, rm22, rm23, rm24, rm25, rm26,
                rm31, rm32, rm33, rm34, rm35, rm36,
                rm41, rm42, rm43, rm44, rm45, rm46,
                rm51, rm52, rm53, rm54, rm55, rm56,
                rm61, rm62, rm63, rm64, rm65, rm66; // 36 RMij parameters
        
        double  t111, t112, t113, t114, t115, t116,
                t121, t122, t123, t124, t125, t126,
                t131, t132, t133, t134, t135, t136,
                t141, t142, t143, t144, t145, t146,
                t151, t152, t153, t154, t155, t156,
                t161, t162, t163, t164, t165, t166,
    
                t211, t212, t213, t214, t215, t216,
                t221, t222, t223, t224, t225, t226,
                t231, t232, t233, t234, t235, t236,
                t241, t242, t243, t244, t245, t246,
                t251, t252, t253, t254, t255, t256,
                t261, t262, t263, t264, t265, t266,
    
                t311, t312, t313, t314, t315, t316,
                t321, t322, t323, t324, t325, t326,
                t331, t332, t333, t334, t335, t336,
                t341, t342, t343, t344, t345, t346,
                t351, t352, t353, t354, t355, t356,
                t361, t362, t363, t364, t365, t366,
    
                t411, t412, t413, t414, t415, t416,
                t421, t422, t423, t424, t425, t426,
                t431, t432, t433, t434, t435, t436,
                t441, t442, t443, t444, t445, t446,
                t451, t452, t453, t454, t455, t456,
                t461, t462, t463, t464, t465, t466,
    
                t511, t512, t513, t514, t515, t516,
                t521, t522, t523, t524, t525, t526,
                t531, t532, t533, t534, t535, t536,
                t541, t542, t543, t544, t545, t546,
                t551, t552, t553, t554, t555, t556,
                t561, t562, t563, t564, t565, t566,
    
                t611, t612, t613, t614, t615, t616,
                t621, t622, t623, t624, t625, t626,
                t631, t632, t633, t634, t635, t636,
                t641, t642, t643, t644, t645, t646,
                t651, t652, t653, t654, t655, t656,
		        t661, t662, t663, t664, t665, t666; // 216 Tijk parameters
	    
    
    };

    template<class BP> // BP is bunch particles?
    struct PropMatrix
    {
        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        const MatrixParams matrix_impl;

        PropMatrix(BP & bp, MatrixParams const& mp)
            : p(bp.parts), masks(bp.masks), mp(mp)
            { }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (masks(i))
            {

                // bend
                FF_algorithm::matrix_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3),
                        dipedge_params.re_2_1, dipedge_params.re_4_3,
                        dipedge_params.te );

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
        const DipedgeParams mp;

        PropDipedgeSimd(BP & bp, DipedgeParams const& mp)
            : p(bp.parts), masks(bp.masks), mp(mp)
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

 
                // bend
                FF_algorithm::dipedge_unit<gsv_t>(
                        p0, p1, p2, p3,
                        dipedge_params.re_2_1,
                        dipedge_params.re_4_3,
                        dipedge_params.te );


                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
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


        FF_algorithm::dipedge_unit(
                x_l, xp_l, y_l, yp_l,
                dipedge_params.re_2_1,
                dipedge_params.re_4_3,
                dipedge_params.te);
        
    }
}

namespace FF_dipedge {

    template <class BunchT>
    inline void
    apply(Lattice_element_slice const& slice, BunchT& bunch)
    {
        using namespace dipedge_impl;

        scoped_simple_timer timer("libFF_dipedge");

        auto const& ele = slice.get_lattice_element();

        MatrixParams mp;

        mp.k1 = ele.get_double_attribute("kick1", 0.0);
        mp.k2 = ele.get_double_attribute("kick2", 0.0);
        mp.k3 = ele.get_double_attribute("kick3", 0.0);
        mp.k4 = ele.get_double_attribute("kick4", 0.0);
        mp.k5 = ele.get_double_attribute("kick5", 0.0);
        mp.k6 = ele.get_double_attribute("kick6", 0.0);

        mp.rm11 = ele.get_double_attribute("rm11", 1.0);
        mp.rm22 = ele.get_double_attribute("rm22", 1.0);
        mp.rm33 = ele.get_double_attribute("rm33", 1.0);
        mp.rm44 = ele.get_double_attribute("rm44", 1.0);
        mp.rm55 = ele.get_double_attribute("rm55", 1.0);
        mp.rm66 = ele.get_double_attribute("rm66", 1.0);

        mp.rm12 = ele.get_double_attribute("rm12", 0.0);
        mp.rm13 = ele.get_double_attribute("rm13", 0.0);
        mp.rm14 = ele.get_double_attribute("rm14", 0.0)
        mp.rm15 = ele.get_double_attribute("rm15", 0.0);
        mp.rm16 = ele.get_double_attribute("rm16", 0.0);
    
        mp.rm21 = ele.get_double_attribute("rm21", 0.0);
        mp.rm23 = ele.get_double_attribute("rm23", 0.0);
        mp.rm24 = ele.get_double_attribute("rm24", 0.0)
        mp.rm25 = ele.get_double_attribute("rm25", 0.0);
        mp.rm26 = ele.get_double_attribute("rm26", 0.0);
    
        mp.rm31 = ele.get_double_attribute("rm31", 0.0);
        mp.rm32 = ele.get_double_attribute("rm32", 0.0);
        mp.rm34 = ele.get_double_attribute("rm34", 0.0)
        mp.rm35 = ele.get_double_attribute("rm35", 0.0);
        mp.rm36 = ele.get_double_attribute("rm36", 0.0);
    
        mp.rm41 = ele.get_double_attribute("rm41", 0.0)
        mp.rm42 = ele.get_double_attribute("rm42", 0.0);
        mp.rm43 = ele.get_double_attribute("rm43", 0.0);
        mp.rm45 = ele.get_double_attribute("rm45", 0.0);
        mp.rm46 = ele.get_double_attribute("rm46", 0.0);
    
        mp.rm51 = ele.get_double_attribute("rm51", 0.0);
        mp.rm52 = ele.get_double_attribute("rm52", 0.0);
        mp.rm53 = ele.get_double_attribute("rm53", 0.0);
        mp.rm54 = ele.get_double_attribute("rm54", 0.0)
        mp.rm56 = ele.get_double_attribute("rm56", 0.0);
    
        mp.rm61 = ele.get_double_attribute("rm61", 0.0);
        mp.rm62 = ele.get_double_attribute("rm62", 0.0);
        mp.rm63 = ele.get_double_attribute("rm63", 0.0);
        mp.rm64 = ele.get_double_attribute("rm64", 0.0)
        mp.rm65 = ele.get_double_attribute("rm65", 0.0);
    
	mp.t111 = ele.get_double_attribute("tm111", 0.0);
	mp.t112 = ele.get_double_attribute("tm112", 0.0);
	mp.t113 = ele.get_double_attribute("tm113", 0.0);
	mp.t114 = ele.get_double_attribute("tm114", 0.0);
	mp.t115 = ele.get_double_attribute("tm115", 0.0);
	mp.t116 = ele.get_double_attribute("tm116", 0.0);

	mp.t121 = ele.get_double_attribute("tm121", 0.0);
	mp.t122 = ele.get_double_attribute("tm122", 0.0);
	mp.t123 = ele.get_double_attribute("tm123", 0.0);
	mp.t124 = ele.get_double_attribute("tm124", 0.0);
	mp.t125 = ele.get_double_attribute("tm125", 0.0);
	mp.t126 = ele.get_double_attribute("tm126", 0.0);

	mp.t131 = ele.get_double_attribute("tm131", 0.0);
	mp.t132 = ele.get_double_attribute("tm132", 0.0);
	mp.t133 = ele.get_double_attribute("tm133", 0.0);
	mp.t134 = ele.get_double_attribute("tm134", 0.0);
	mp.t135 = ele.get_double_attribute("tm135", 0.0);
	mp.t136 = ele.get_double_attribute("tm136", 0.0);

	mp.t131 = ele.get_double_attribute("tm131", 0.0);
	mp.t132 = ele.get_double_attribute("tm132", 0.0);
	mp.t133 = ele.get_double_attribute("tm133", 0.0);
	mp.t134 = ele.get_double_attribute("tm134", 0.0);
	mp.t135 = ele.get_double_attribute("tm135", 0.0);
	mp.t136 = ele.get_double_attribute("tm136", 0.0);

	dipedge_params.re_2_1 = h * tanedg;
        dipedge_params.re_4_3 = -h * tan(psip);


        Reference_particle& ref_l = bunch.get_design_reference_particle();
        Reference_particle const& ref_b = bunch.get_reference_particle();

        dipedge_params.scale =
            ref_l.get_momentum() /
            (ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]));

        // lattice reference
        double pref_l = ref_l.get_momentum();
        double brho_l = pref_l / PH_CNV_brho_to_p;
        double m_l = ref_l.get_mass();
        int charge_l = ref_l.get_charge();

        // the particle dpop is with respect to this momentum which goes
        // with the bunch
        double pref_b = ref_b.get_momentum();
        double brho_b = pref_b / PH_CNV_brho_to_p;
        double m_b = bunch.get_mass();

        // common
        dipedge_params.pref_b = pref_b;
        dipedge_params.m_b = m_b;

        dipedge_params.ce1 = cos(-edge);
        dipedge_params.se1 = sin(-edge);

        using namespace Kokkos;
        using exec = typename BunchT::exec_space;

        auto apply = [&](ParticleGroup pg) {
            auto bp = bunch.get_bunch_particles(pg);
            if (!bp.num_valid()) return;

            // dipedge is a 0 length element with no multipole moments
            // (currently)
#if LIBFF_USE_GSV
            prop_reference(ref_l, dipedge_params);

            auto range = RangePolicy<exec>(0, bp.size_in_gsv());
            PropDipedgeSimd<typename BunchT::bp_t> dipedge(bp, dipedge_params);
            Kokkos::parallel_for(range, dipedge);

#else  // LIBFF_USE_GSV
            auto range = RangePolicy<exec>(0, bp.size());
            PropDipedge<typename BunchT::bp_t> dipedge(bp, dipedge_params);
            Kokkos::parallel_for(range, dipedge);
#endif // LIBFF_USE_GSV
        };

        apply(ParticleGroup::regular);
        apply(ParticleGroup::spectator);
        // dipedge has 0 length so there is no need to increment trajectory here
        // bunch.get_reference_particle().increment_trajectory(0);
        Kokkos::fence();
    }
}

#endif // FF_MATRIX_H


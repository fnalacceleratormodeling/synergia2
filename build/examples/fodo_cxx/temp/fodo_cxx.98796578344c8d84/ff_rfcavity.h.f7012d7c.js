var data = {lines:[
{"lineNum":"    1","line":"#ifndef FF_RFCAVITY_H"},
{"lineNum":"    2","line":"#define FF_RFCAVITY_H"},
{"lineNum":"    3","line":""},
{"lineNum":"    4","line":"#include \"synergia/foundation/math_constants.h\""},
{"lineNum":"    5","line":"#include \"synergia/libFF/ff_algorithm.h\""},
{"lineNum":"    6","line":"#include \"synergia/libFF/ff_patterned_propagator.h\""},
{"lineNum":"    7","line":"#include \"synergia/utils/simple_timer.h\""},
{"lineNum":"    8","line":""},
{"lineNum":"    9","line":""},
{"lineNum":"   10","line":"namespace rfcavity_impl"},
{"lineNum":"   11","line":"{"},
{"lineNum":"   12","line":"    struct RFCavityParams"},
{"lineNum":"   13","line":"    {"},
{"lineNum":"   14","line":"        int nh;"},
{"lineNum":"   15","line":"        double mhp[12];"},
{"lineNum":"   16","line":""},
{"lineNum":"   17","line":"        double length;"},
{"lineNum":"   18","line":"        double str;"},
{"lineNum":"   19","line":"        double phi_s;"},
{"lineNum":"   20","line":"        double w_rf;"},
{"lineNum":"   21","line":""},
{"lineNum":"   22","line":"        double pref_b;"},
{"lineNum":"   23","line":"        double m_b;"},
{"lineNum":"   24","line":"        double new_pref_b;"},
{"lineNum":"   25","line":""},
{"lineNum":"   26","line":"        double ref_cdt_1;"},
{"lineNum":"   27","line":"        double ref_cdt_2;"},
{"lineNum":"   28","line":"    };"},
{"lineNum":"   29","line":""},
{"lineNum":"   30","line":"    template<class BunchT>"},
{"lineNum":"   31","line":"    struct PropThinRFCavity","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"   32","line":"    {"},
{"lineNum":"   33","line":"        using gsv_t = typename BunchT::gsv_t;"},
{"lineNum":"   34","line":""},
{"lineNum":"   35","line":"        typename BunchT::bp_t::parts_t p;"},
{"lineNum":"   36","line":"        typename BunchT::bp_t::const_masks_t m;"},
{"lineNum":"   37","line":"        const RFCavityParams rp;"},
{"lineNum":"   38","line":""},
{"lineNum":"   39","line":"        KOKKOS_INLINE_FUNCTION"},
{"lineNum":"   40","line":"        void operator()(const int idx) const"},
{"lineNum":"   41","line":"        {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   42","line":"            int i = idx * gsv_t::size();"},
{"lineNum":"   43","line":"            int mask = 0;"},
{"lineNum":"   44","line":"            for(int x=i; x<i+gsv_t::size(); ++x) mask |= m(x);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   45","line":""},
{"lineNum":"   46","line":"            if (mask)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   47","line":"            {"},
{"lineNum":"   48","line":"                gsv_t p1(&p(i, 1));"},
{"lineNum":"   49","line":"                gsv_t p3(&p(i, 3));"},
{"lineNum":"   50","line":"                gsv_t p4(&p(i, 4));"},
{"lineNum":"   51","line":"                gsv_t p5(&p(i, 5));"},
{"lineNum":"   52","line":""},
{"lineNum":"   53","line":"                FF_algorithm::thin_rfcavity_unit(","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   54","line":"                        p1, p3, p4, p5,"},
{"lineNum":"   55","line":"                        rp.w_rf, rp.str, rp.phi_s, rp.m_b, rp.pref_b, rp.new_pref_b,","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   56","line":"                        rp.mhp, rp.nh);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   57","line":""},
{"lineNum":"   58","line":"                p1.store(&p(i, 1));"},
{"lineNum":"   59","line":"                p3.store(&p(i, 3));"},
{"lineNum":"   60","line":"                p4.store(&p(i, 4));"},
{"lineNum":"   61","line":"                p5.store(&p(i, 5));"},
{"lineNum":"   62","line":"            }"},
{"lineNum":"   63","line":"        }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   64","line":"    };"},
{"lineNum":"   65","line":""},
{"lineNum":"   66","line":""},
{"lineNum":"   67","line":"    template<class BunchT>"},
{"lineNum":"   68","line":"    struct PropRFCavity","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"   69","line":"    {"},
{"lineNum":"   70","line":"        using gsv_t = typename BunchT::gsv_t;"},
{"lineNum":"   71","line":""},
{"lineNum":"   72","line":"        typename BunchT::bp_t::parts_t p;"},
{"lineNum":"   73","line":"        typename BunchT::bp_t::const_masks_t m;"},
{"lineNum":"   74","line":"        const RFCavityParams rp;"},
{"lineNum":"   75","line":""},
{"lineNum":"   76","line":"        KOKKOS_INLINE_FUNCTION"},
{"lineNum":"   77","line":"        void operator()(const int idx) const"},
{"lineNum":"   78","line":"        {","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   79","line":"            int i = idx * gsv_t::size();","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   80","line":"            int mask = 0;"},
{"lineNum":"   81","line":"            for(int x=i; x<i+gsv_t::size(); ++x) mask |= m(x);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   82","line":""},
{"lineNum":"   83","line":"            if (mask)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   84","line":"            {"},
{"lineNum":"   85","line":"                gsv_t p0(&p(i, 0));"},
{"lineNum":"   86","line":"                gsv_t p1(&p(i, 1));"},
{"lineNum":"   87","line":"                gsv_t p2(&p(i, 2));"},
{"lineNum":"   88","line":"                gsv_t p3(&p(i, 3));"},
{"lineNum":"   89","line":"                gsv_t p4(&p(i, 4));"},
{"lineNum":"   90","line":"                gsv_t p5(&p(i, 5));"},
{"lineNum":"   91","line":""},
{"lineNum":"   92","line":"                FF_algorithm::drift_unit(","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   93","line":"                        p0, p1, p2, p3, p4, p5,"},
{"lineNum":"   94","line":"                        0.5 * rp.length, rp.pref_b, rp.m_b, rp.ref_cdt_1);","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"   95","line":""},
{"lineNum":"   96","line":"                FF_algorithm::thin_rfcavity_unit(","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   97","line":"                        p1, p3, p4, p5,"},
{"lineNum":"   98","line":"                        rp.w_rf, rp.str, rp.phi_s, rp.m_b, rp.pref_b, rp.new_pref_b,","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   99","line":"                        rp.mhp, rp.nh);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  100","line":""},
{"lineNum":"  101","line":"                FF_algorithm::drift_unit(","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  102","line":"                        p0, p1, p2, p3, p4, p5,"},
{"lineNum":"  103","line":"                        0.5 * rp.length, rp.pref_b, rp.m_b, rp.ref_cdt_2);","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  104","line":""},
{"lineNum":"  105","line":"                p0.store(&p(i, 0));"},
{"lineNum":"  106","line":"                p1.store(&p(i, 1));"},
{"lineNum":"  107","line":"                p2.store(&p(i, 2));"},
{"lineNum":"  108","line":"                p3.store(&p(i, 3));"},
{"lineNum":"  109","line":"                p4.store(&p(i, 4));"},
{"lineNum":"  110","line":"                p5.store(&p(i, 5));"},
{"lineNum":"  111","line":"            }"},
{"lineNum":"  112","line":"        }","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  113","line":"    };"},
{"lineNum":"  114","line":""},
{"lineNum":"  115","line":""},
{"lineNum":"  116","line":"}"},
{"lineNum":"  117","line":""},
{"lineNum":"  118","line":"namespace FF_rfcavity"},
{"lineNum":"  119","line":"{"},
{"lineNum":"  120","line":"    template<class BunchT>"},
{"lineNum":"  121","line":"    void apply(Lattice_element_slice const& slice, BunchT& bunch)"},
{"lineNum":"  122","line":"    {","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  123","line":"        using namespace rfcavity_impl;"},
{"lineNum":"  124","line":""},
{"lineNum":"  125","line":"        scoped_simple_timer timer(\"libFF_rfcavity\");"},
{"lineNum":"  126","line":""},
{"lineNum":"  127","line":"        RFCavityParams rp;"},
{"lineNum":"  128","line":""},
{"lineNum":"  129","line":"        rp.length = slice.get_right() - slice.get_left();","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  130","line":"        Lattice_element const & elm = slice.get_lattice_element();"},
{"lineNum":"  131","line":""},
{"lineNum":"  132","line":"        int    harmonic_number = elm.get_double_attribute(\"harmon\", -1.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  133","line":"        double            volt = elm.get_double_attribute(\"volt\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  134","line":"        double             lag = elm.get_double_attribute(\"lag\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  135","line":"        double           shunt = elm.get_double_attribute(\"shunt\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  136","line":"        // freq is the synchronous frequency of the cavity in MHz"},
{"lineNum":"  137","line":"        double            freq = elm.get_double_attribute(\"freq\", -1.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  138","line":"        // delta_freq is the frequency offset from synchronous in MHz like freq"},
{"lineNum":"  139","line":"        double      delta_freq = elm.get_double_attribute(\"delta_freq\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  140","line":""},
{"lineNum":"  141","line":"        // harmonics,"},
{"lineNum":"  142","line":"        // mhp[h*3]   = harmonic multiple"},
{"lineNum":"  143","line":"        // mhp[h*3+1] = relative strength"},
{"lineNum":"  144","line":"        // mhp[h*3+2] = phase shift"},
{"lineNum":"  145","line":"        // nh = number of harmonics"},
{"lineNum":"  146","line":"        rp.nh = 1;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  147","line":"        rp.mhp[0] = 1.0;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  148","line":"        rp.mhp[1] = 1.0;"},
{"lineNum":"  149","line":"        rp.mhp[2] = 0.0;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  150","line":""},
{"lineNum":"  151","line":"        // harmonic multiple = 2"},
{"lineNum":"  152","line":"        if (elm.has_double_attribute(\"h2_str\"))","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  153","line":"        {"},
{"lineNum":"  154","line":"            rp.mhp[rp.nh*3+0] = 2;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  155","line":"            rp.mhp[rp.nh*3+1] = elm.get_double_attribute(\"h2_str\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  156","line":"            rp.mhp[rp.nh*3+2] = elm.get_double_attribute(\"h2_phase\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  157","line":"            ++rp.nh;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  158","line":"        }"},
{"lineNum":"  159","line":""},
{"lineNum":"  160","line":"        // harmonic multiple = 3"},
{"lineNum":"  161","line":"        if (elm.has_double_attribute(\"h3_str\"))","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  162","line":"        {"},
{"lineNum":"  163","line":"            rp.mhp[rp.nh*3+0] = 3;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  164","line":"            rp.mhp[rp.nh*3+1] = elm.get_double_attribute(\"h3_str\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  165","line":"            rp.mhp[rp.nh*3+2] = elm.get_double_attribute(\"h3_phase\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  166","line":"            ++rp.nh;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  167","line":"        }"},
{"lineNum":"  168","line":""},
{"lineNum":"  169","line":"        // harmonic multiple = 4"},
{"lineNum":"  170","line":"        if (elm.has_double_attribute(\"h4_str\"))","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  171","line":"        {"},
{"lineNum":"  172","line":"            rp.mhp[rp.nh*3+0] = 4;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  173","line":"            rp.mhp[rp.nh*3+1] = elm.get_double_attribute(\"h4_str\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  174","line":"            rp.mhp[rp.nh*3+2] = elm.get_double_attribute(\"h4_phase\", 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  175","line":"            ++rp.nh;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  176","line":"        }"},
{"lineNum":"  177","line":""},
{"lineNum":"  178","line":"        rp.str = volt * 1.0e-3;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  179","line":""},
{"lineNum":"  180","line":"        // keep lag within the range of [0, 1)."},
{"lineNum":"  181","line":"        while (lag < 0.0)  { lag += 1.0; }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  182","line":"        while (lag >= 1.0) { lag -= 1.0; }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  183","line":""},
{"lineNum":"  184","line":"        //elm.set_double_attribute(\"lag\", lag);"},
{"lineNum":"  185","line":"        rp.phi_s = 2.0 * mconstants::pi * lag;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  186","line":"        rp.w_rf  = 2.0 * mconstants::pi * (freq + delta_freq) * 1.0e6;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  187","line":""},
{"lineNum":"  188","line":"        Reference_particle & ref_l = bunch.get_design_reference_particle();"},
{"lineNum":"  189","line":"        Reference_particle & ref_b = bunch.get_reference_particle();"},
{"lineNum":"  190","line":""},
{"lineNum":"  191","line":"        // The bunch particles momentum is with respect to the bunch reference particle"},
{"lineNum":"  192","line":"        rp.pref_b = ref_b.get_momentum();","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  193","line":"        rp.m_b = bunch.get_mass();","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  194","line":""},
{"lineNum":"  195","line":"        rp.new_pref_b = FF_algorithm::thin_rfcavity_pnew(","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  196","line":"                rp.pref_b, rp.m_b, rp.str, rp.phi_s);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  197","line":""},
{"lineNum":"  198","line":"        // reference_cdt uses the lattice reference particle"},
{"lineNum":"  199","line":"        // double reference_cdt = get_reference_cdt(length, ref_l);"},
{"lineNum":"  200","line":"        double ref_l_x    = ref_l.get_state()[Bunch::x];","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  201","line":"        double ref_l_xp   = ref_l.get_state()[Bunch::xp];","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  202","line":"        double ref_l_y    = ref_l.get_state()[Bunch::y];","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  203","line":"        double ref_l_yp   = ref_l.get_state()[Bunch::yp];","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  204","line":"        double ref_l_cdt  = 0.0;"},
{"lineNum":"  205","line":"        double ref_l_dpop = ref_l.get_state()[Bunch::dpop];","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  206","line":""},
{"lineNum":"  207","line":"        double ref_l_p = ref_l.get_momentum();","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  208","line":"        double ref_l_m = ref_l.get_mass();","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  209","line":""},
{"lineNum":"  210","line":"        // double new_ref_l_p ="},
{"lineNum":"  211","line":"        //     FF_algorithm::thin_rfcavity_pnew(ref_l_p, ref_l_m, str, phi_s);"},
{"lineNum":"  212","line":""},
{"lineNum":"  213","line":"        // first half drift"},
{"lineNum":"  214","line":"        FF_algorithm::drift_unit("},
{"lineNum":"  215","line":"                ref_l_x, ref_l_xp, ref_l_y, ref_l_yp, ref_l_cdt, ref_l_dpop,"},
{"lineNum":"  216","line":"                0.5 * rp.length, ref_l_p, ref_l_m, 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  217","line":""},
{"lineNum":"  218","line":"        double total_ref_cdt = ref_l_cdt;"},
{"lineNum":"  219","line":"        rp.ref_cdt_1 = ref_l_cdt;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  220","line":"        ref_l_cdt = 0.0;"},
{"lineNum":"  221","line":""},
{"lineNum":"  222","line":"        // do not give it the new_ref_l_p because the xp and yp dont get scaled, the momentum"},
{"lineNum":"  223","line":"        // of the lattice reference particle remains unchanged, only the dpop of the state"},
{"lineNum":"  224","line":"        // has been changed"},
{"lineNum":"  225","line":"        FF_algorithm::thin_rfcavity_unit("},
{"lineNum":"  226","line":"                ref_l_xp, ref_l_yp, ref_l_cdt, ref_l_dpop,"},
{"lineNum":"  227","line":"                rp.w_rf, rp.str, rp.phi_s, ref_l_m, ref_l_p, ref_l_p,","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  228","line":"                rp.mhp, rp.nh );"},
{"lineNum":"  229","line":""},
{"lineNum":"  230","line":"        // second half drift"},
{"lineNum":"  231","line":"        FF_algorithm::drift_unit("},
{"lineNum":"  232","line":"                ref_l_x, ref_l_xp, ref_l_y, ref_l_yp, ref_l_cdt, ref_l_dpop,"},
{"lineNum":"  233","line":"                0.5 * rp.length, ref_l_p, ref_l_m, 0.0);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  234","line":""},
{"lineNum":"  235","line":"        total_ref_cdt += ref_l_cdt;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  236","line":"        rp.ref_cdt_2 = ref_l_cdt;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  237","line":""},
{"lineNum":"  238","line":"        // save the state"},
{"lineNum":"  239","line":"        ref_l.set_state(ref_l_x, ref_l_xp, ref_l_y, ref_l_yp, total_ref_cdt, ref_l_dpop);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  240","line":""},
{"lineNum":"  241","line":"        // bunch particles"},
{"lineNum":"  242","line":"        auto apply = [&](ParticleGroup pg) {","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  243","line":"            if (!bunch.get_local_num(pg)) return;","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  244","line":""},
{"lineNum":"  245","line":"            auto parts = bunch.get_local_particles(pg);"},
{"lineNum":"  246","line":"            auto masks = bunch.get_local_particle_masks(pg);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  247","line":""},
{"lineNum":"  248","line":"            using exec = typename BunchT::exec_space;"},
{"lineNum":"  249","line":"            auto range = Kokkos::RangePolicy<exec>(0, bunch.size_in_gsv(pg));","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  250","line":""},
{"lineNum":"  251","line":"            if (close_to_zero(rp.length))","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  252","line":"            {"},
{"lineNum":"  253","line":"                PropThinRFCavity<BunchT> rfcavity{parts, masks, rp};","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  254","line":"                Kokkos::parallel_for(range, rfcavity);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  255","line":"            }"},
{"lineNum":"  256","line":"            else"},
{"lineNum":"  257","line":"            {"},
{"lineNum":"  258","line":"                PropRFCavity<BunchT> rfcavity{parts, masks, rp};","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  259","line":"                Kokkos::parallel_for(range, rfcavity);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  260","line":"            }"},
{"lineNum":"  261","line":"        };","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  262","line":""},
{"lineNum":"  263","line":"        apply(ParticleGroup::regular);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  264","line":"        apply(ParticleGroup::spectator);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  265","line":""},
{"lineNum":"  266","line":"        // updated four momentum"},
{"lineNum":"  267","line":"        Four_momentum fm = ref_b.get_four_momentum();","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  268","line":"        fm.set_momentum(rp.new_pref_b);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  269","line":""},
{"lineNum":"  270","line":"        // update the bunch reference particle with the updated ref_p"},
{"lineNum":"  271","line":"        ref_b.set_four_momentum(fm);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  272","line":"        ref_b.increment_trajectory(rp.length);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  273","line":"    }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  274","line":"}"},
{"lineNum":"  275","line":""},
{"lineNum":"  276","line":"#endif // FF_RFCAVITY_H"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 87, "covered" : 0,};
var merged_data = [];
